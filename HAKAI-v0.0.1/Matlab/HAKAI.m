function res = HAKAI(varargin)

% HAKAI - A 3-dimensional finite element program
% Copyright (c) 2024 Yozo Yugen
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see http://www.gnu.org/licenses/.


    fname = sprintf("matlab_bug_report.txt");
    bug_report = fopen(fname,"w");


    varargin
       
    [ MODEL] = readInpFile(varargin{1});
    nNode = MODEL.nNode;
    coordmat = MODEL.coordmat;
    nElement = MODEL.nElement;

    elementmat = MODEL.elementmat;
    element_material = MODEL.element_material;
    element_instance = MODEL.element_instance;
    contact_flag = MODEL.contact_flag
    
    flag_fracture = 0;
    
    mass_scaling = MODEL.mass_scaling;
    d_time = MODEL.d_time * sqrt(mass_scaling);
    end_time = MODEL.end_time;
    time_num = floor(end_time / d_time)
    C = 0.0;


    force_node = [];
    force_dof = [];
    force_v = [];        
   
    
  
    close all
    drawElement(coordmat, elementmat, 1);


    %--- Material property 
    for i = 1 : length(MODEL.MATERIAL)
        young = MODEL.MATERIAL(i).young;
        poisson = MODEL.MATERIAL(i).poisson;
        G = young / 2 / (1+poisson);
        MODEL.MATERIAL(i).G = G;
    
        %Dmat = zeros(6,6);
        d1 = (1.0-poisson);
        d2 = poisson;
        d3 = (1.0-2.0*poisson)/2.0;
        Dmat = young / (1.0+poisson) / (1-2.0*poisson)*...
           [ d1 d2 d2 0  0  0
             d2 d1 d2 0  0  0
             d2 d2 d1 0  0  0
             0  0  0  d3 0  0
             0  0  0  0  d3 0
             0  0  0  0  0  d3];
        MODEL.MATERIAL(i).Dmat = Dmat;

        if length( MODEL.MATERIAL(i).failure_stress ) > 0
            flag_fracture = 1;
            MODEL.MATERIAL(i).failure_stress
        end

        if size( MODEL.MATERIAL(i).ductile, 1) > 0
            flag_fracture = 1;        
        end
        
    end

    
    
    %--- set variable
    fn = nNode * 3;
    integ_num = 8;
    %integ_num = 1;

    %--- shape func
    Pusai_mat = zeros(3,8,integ_num);
    Pusai_mat = cal_Pusai_hexa(integ_num);
    
    % elementdofmat = zeros(8*3,nElement);
    % for e = 1 : nElement
    %     for j = 1 : 8
    %         elementdofmat(j*3-2,e) = elementmat(j,e)*3-2;
    %         elementdofmat(j*3-1,e) = elementmat(j,e)*3-1;
    %         elementdofmat(j*3,e) = elementmat(j,e)*3;
    %     end
    % end

    elementVolume = zeros(nElement,1);
    for e = 1 : nElement
        e_position = coordmat(:, elementmat(:,e));
        V = 0;
        for i = 1 : integ_num
            J = Pusai_mat(:,:,i) * e_position';
            V = V + my3det(J);
        end
        elementVolume(e) = V;
    end
    elementVolume.';

    diag_M = zeros(fn, 1);
    for i = 1 : nElement
        density = MODEL.MATERIAL( element_material(i) ).density;
        node_mass = density * elementVolume(i) / 8;
        diag_M( (elementmat(:,i)-1)*3 + 1 ) = diag_M( (elementmat(:,i)-1)*3 + 1 ) + node_mass;
        diag_M( (elementmat(:,i)-1)*3 + 2 ) = diag_M( (elementmat(:,i)-1)*3 + 2 ) + node_mass;
        diag_M( (elementmat(:,i)-1)*3 + 3 ) = diag_M( (elementmat(:,i)-1)*3 + 3 ) + node_mass;
    end

    diag_M = diag_M * mass_scaling;
    diag_C = diag_M * C;

    % for i = 1 : fn
    %     fprintf(bug_report, "%d, %.16e\n", i, diag_M(i) );
    % end


    position = coordmat;
    disp = zeros(fn, 1);
    disp_pre = zeros(fn, 1);
    velo = zeros(fn, 1);

    %--- Initial Velocity ---%
    for i = 1 : length(MODEL.IC)
        disp_pre( MODEL.IC(i).dof ) = -MODEL.IC(i).value * d_time * 1;
        %disp_pre.'
        velo( MODEL.IC(i).dof ) = MODEL.IC(i).value;
    end

    %--- Contact ---%
    %c_pair = [];
    c_INSTANCE = [];
    % c_INSTANCE(1:length(MODEL.INSTANCE)) =...
    %     struct( 'part_id', 0, 'material_id', 0, ...
    %             'nNode', 0, 'node_offset', 0, 'element_offset', 0, ...
    %             'c_triangles', [],...
    %             'c_triangles_eleid', [],...
    %             'c_nodes', [] );

    if contact_flag >= 1
        tic;

        for i = 1 : length(MODEL.INSTANCE)  
            [faces, faces_eleid, sorted_faces] = get_element_face(MODEL, i);
            MODEL.INSTANCE(i).surfaces = faces;
            MODEL.INSTANCE(i).surfaces_eleid = faces_eleid;
            MODEL.INSTANCE(i).sorted_surfaces = sorted_faces;
        end


        if length(MODEL.CP) == 0  % ALL EXTERIOR

            ni = length(MODEL.INSTANCE);
            
            if ni > 1
                c = 1;
                for i = 1 : ni
                    js = i+1;
                    %js = i;  % -> include self-contact
                    if contact_flag == 2
                        js = i;
                    end

                    for j = js : ni
                        MODEL.CP(c).instance_id_i = i;
                        MODEL.CP(c).instance_id_j = j;
                        MODEL.CP(c).elements_i = 1 : MODEL.INSTANCE(i).nElement;
                        MODEL.CP(c).elements_j = 1 : MODEL.INSTANCE(j).nElement;
                        c = c+1;
                    end
                end    

            else % ni == 1
                MODEL.CP(1).instance_id_i = 1;
                MODEL.CP(1).instance_id_j = 1;
                MODEL.CP(1).elements_i = 1 : MODEL.INSTANCE(1).nElement;
                MODEL.CP(1).elements_j = 1 : MODEL.INSTANCE(1).nElement;
            end

        else

        end
        %c_pair
    
        for i = 1 : length(MODEL.CP)  

            array_element = 1 : MODEL.INSTANCE(MODEL.CP(i).instance_id_i).nElement;
            [c_triangles, c_triangles_eleid, c_nodes] = ...
                get_surface_triangle(MODEL.INSTANCE(MODEL.CP(i).instance_id_i), array_element, MODEL.CP(i).elements_i);
            MODEL.CP(i).c_triangles_i = c_triangles;
            MODEL.CP(i).c_triangles_eleid_i = c_triangles_eleid;
            MODEL.CP(i).c_nodes_i = c_nodes;
            
            array_element = 1 : MODEL.INSTANCE(MODEL.CP(i).instance_id_j).nElement;
            [c_triangles, c_triangles_eleid, c_nodes] = ...
                get_surface_triangle(MODEL.INSTANCE(MODEL.CP(i).instance_id_j), array_element, MODEL.CP(i).elements_j);
            MODEL.CP(i).c_triangles_j = c_triangles;
            MODEL.CP(i).c_triangles_eleid_j = c_triangles_eleid;
            MODEL.CP(i).c_nodes_j = c_nodes;

            MODEL.CP(i)
        end

        
        'Contact set'
        toc

        % tic
        %     for i = 1 : nNode-1
        %         for j = i+1 : nNode
        %             %v = vecnorm(coordmat(:,i) - coordmat(:,j)); %-> same
        %             v = norm(coordmat(:,i) - coordmat(:,j));
        %         end
        %     end
        %     'Neighbor nodes'
        % toc

        % tic
        %     nodeNormal = zeros(3, nNode);
        %     elementNormal = zeros(3, nElement);
        %     for i = 1 : nElement
        %         nd = elementmat(:,i);
        %         ii = nd([1 2 3 4 1]);
        %         line( coordmat(1,ii),  coordmat(2,ii), coordmat(3,ii), 'color', 'b' )
        % 
        %         v1 = coordmat(ii(2), :) - coordmat(ii(1), :);
        %         v4 = coordmat(ii(4), :) - coordmat(ii(1), :);
        %         n = cross(v1,v4);
        %         n = n / vecnorm(n);
        %         elementNormal(:,i) = n;
        % 
        %         nodeNormal(ii(1), :) = nodeNormal(ii(1), :) + n;
        %         nodeNormal(ii(2), :) = nodeNormal(ii(2), :) + n;
        %         nodeNormal(ii(3), :) = nodeNormal(ii(3), :) + n;
        %         nodeNormal(ii(4), :) = nodeNormal(ii(4), :) + n;
        %     end
        %     'node normal'
        % toc

    end


    elementSize = zeros(nElement,3);
    for e = 1 : nElement
        e_pos = coordmat(:, elementmat(:,e));
        L1 = norm( e_pos(:,1) - e_pos(:,2) );
        L2 = norm( e_pos(:,1) - e_pos(:,4) );
        L3 = norm( e_pos(:,1) - e_pos(:,5) );
        elementSize(e,:) = [L1 L2 L3];
    end
    elementMinSize = min(min(elementSize))
    elementMaxSize = max(max(elementSize))
    d_max = 0;
    d_node = zeros(nNode,1);
    d_disp_norm = zeros(nNode,1);
    c_force3 = zeros(3, nNode);

    %--- Variable ---%
    Q = zeros(fn, 1);
    % integ_strain = zeros( nElement * integ_num, 6); % row-base is slower
    integ_stress = zeros( 6, nElement * integ_num);
    integ_strain = zeros( 6, nElement * integ_num);
    integ_plastic_strain = zeros( 6, nElement * integ_num);
    integ_eq_plastic_strain = zeros( nElement * integ_num, 1);
    integ_triax_stress = zeros( nElement * integ_num, 1);
    element_flag = ones(nElement, 1);
    integ_flag = ones(nElement, integ_num);
  
    integ_yield_stress = zeros( nElement * integ_num, 1);
    for i = 1 : nElement
        pp = MODEL.MATERIAL(element_material(i)).plastic;
        if length( pp ) >0
            integ_yield_stress((1:integ_num)+(i-1)*integ_num) = pp(1,1);
        end
    end

    element_ctr = zeros(3, nElement);
    element_ctr = cal_element_ctr(coordmat, elementmat);


    output_num = 100
    d_out = floor(time_num / output_num)
    output_data(1:output_num) = struct('disp',zeros(fn, 1),...
                               'integ_stress',zeros(6, nElement * integ_num),...  
                               'integ_strain',zeros(6, nElement * integ_num),...  
                               'integ_plastic_strain', zeros(6, nElement * integ_num),... 
                               'integ_eq_plastic_strain', zeros( nElement * integ_num, 1),...
                               'integ_triax_stress', zeros( nElement * integ_num, 1),...
                               'element_flag',ones(nElement, 1) );

    [node_value] = cal_node_stress_strain(nNode, elementmat, integ_num, output_data(1));
    write_vtk(0, coordmat, elementmat, output_data(1).element_flag, output_data(1).disp, velo, node_value);

    i_out = 1;
    clstr_num = 0;

    output_ = zeros( time_num,1);

    tic;
    for t = 1 : time_num

      if rem(t,100) == 0
        % t * d_time
        str = sprintf('time = %.4e / %.4e\n', t*d_time, end_time);
        clstr = repmat('\b',1, clstr_num);
        fprintf([clstr str])
        clstr_num = length(str);
      end

      %set force
      external_force = zeros(fn, 1);
      external_force( force_dof ) = force_v;
      
      if contact_flag >= 1
          %fprintf('contact\n')
          %tic
          [c_force3, d_node] = cal_contact_force(MODEL.CP, position, velo, diag_M, elementMinSize, elementMaxSize, d_max*1, d_node, ...
              MODEL.INSTANCE, MODEL.MATERIAL, element_flag, elementmat, bug_report, t*d_time); % -> even
          %toc
          external_force = external_force + reshape(c_force3, fn,1);
      end


        % update position
        disp_new = 1.0 ./ (diag_M/d_time^2.0 + diag_C/2.0/d_time) .* ( external_force - Q + diag_M/d_time^2.0 .* (2*disp - disp_pre) + diag_C/2.0/d_time.*disp_pre );
        
        % Boundary Conditions
        %tic;
        for i = 1 : length(MODEL.BC)
            amp = 1;
            if length(MODEL.BC(i).amplitude) > 0
                time_index = 1;
                current_time = t * d_time;
                a_t = MODEL.BC(i).amplitude.time;
                a_v = MODEL.BC(i).amplitude.value;
                for j = 1 : length(a_t) - 1
                    if current_time >= a_t(j) && ...
                            current_time <= a_t(j+1)
                        time_index = j;
                        break
                    end
                end
                amp = a_v(time_index) + (a_v(time_index+1) - a_v(time_index)) ...
                      * (current_time - a_t(time_index))/(a_t(time_index+1)-a_t(time_index));
            end

            for j = 1 : length(MODEL.BC(i).dof)
                dof = MODEL.BC(i).dof{j};
                v = MODEL.BC(i).value{j};
                %disp_new( dof ) = v * amp;
                for k = 1 : length(dof)
                    disp_new( dof(k) ) = v * amp;
                end
            end
        end
        %toc

        d_disp = disp_new - disp;
        disp_pre = disp;
        disp = disp_new;
        velo = d_disp / d_time;

        d_disp3 = reshape(d_disp,3,nNode);
        for i = 1 : nNode
            d_disp_norm(i) = my3norm(d_disp3(:,i));
        end
        d_max = max(d_disp_norm);
        %fprintf('d_max: %.10e\n', d_max)
        %output_ = [output_ d_max];
        output_(t) = d_max;
        

        position = coordmat + reshape(disp,3,nNode);

        element_ctr = cal_element_ctr(position, elementmat);


      %--- cal element stress, strain,     
      %                                    
      [d_integ_stress, d_integ_strain, d_integ_yield_stress, d_integ_plastic_strain, d_integ_eq_plastic_strain, Q] = ...
          cal_stress_hexa(position, d_disp, elementmat, element_flag, integ_num, Pusai_mat, integ_stress, ...
                          MODEL.MATERIAL, element_material, integ_yield_stress, integ_eq_plastic_strain, elementMinSize);

      %--- Updated integ values            
      integ_stress = integ_stress + d_integ_stress;
      integ_strain = integ_strain + d_integ_strain;
      integ_plastic_strain = integ_plastic_strain + d_integ_plastic_strain;
      integ_eq_plastic_strain = integ_eq_plastic_strain + d_integ_eq_plastic_strain;
      integ_yield_stress = integ_yield_stress + d_integ_yield_stress;
      integ_triax_stress = cal_triax_stress( integ_stress );

      % for i = 1 : nElement * integ_num
      %   integ_triax_stress(i) = cal_triax_stress_e( integ_stress(:,i) ); %-> even
      % end


      deleted_element = [];
      if flag_fracture == 1

          %--- fracture stress ---%
          % for e = 1 : nElement*integ_num
          %   e_index = fix((e-0.1)/integ_num) + 1;
          %   mat_id = element_material(e_index);
          %   fracture_stress = MODEL.MATERIAL(mat_id).failure;
          %   if length( fracture_stress ) > 0
          %       if integ_yield_stress(e) > fracture_stress && element_flag(e_index) == 1
          %           element_flag( e_index ) = 0;
          %           integ_stress( [1:integ_num]+(e_index-1)*integ_num, :) = 0;
          %           integ_strain( [1:integ_num]+(e_index-1)*integ_num, :) = 0;
          %       end
          %   end
          % end

          %--- fracture strain ---%
          for i = 1 : nElement
            mat_id = element_material(i);
            
            nd =size(MODEL.MATERIAL(mat_id).ductile, 1);
            if nd > 0

                v_e = sum( integ_eq_plastic_strain((1:integ_num)+(i-1)*integ_num) ) / integ_num;
                t_e = sum( integ_triax_stress((1:integ_num)+(i-1)*integ_num) ) / integ_num;
                if t_e < 0
                    continue
                end

                ductile_ = MODEL.MATERIAL(mat_id).ductile;
                fr_e = ductile_(nd,1);
                for j = 1 : nd-1
                    if t_e >= ductile_(j,2) && t_e < ductile_(j+1, 2)
                        fr_e = ductile_(j,1) + (ductile_(j+1,1) - ductile_(j,1)) / (ductile_(j+1,2) - ductile_(j,2)) * (t_e - ductile_(j,2));
                        break
                    end
                end


                if v_e >= fr_e  &&  element_flag(i) == 1
                    % i
                    % v_e
                    % t_e
                    % fr_e
                    element_flag( i ) = 0;
                    integ_stress(:, (1:integ_num)+(i-1)*integ_num) = 0;
                    integ_strain(:, (1:integ_num)+(i-1)*integ_num) = 0;
                    deleted_element = [deleted_element i];
                    %[sum(element_flag)  nElement]
                    str = sprintf('Elememts: %d / %d\n', sum(element_flag),  nElement);
                    fprintf(str)
                    clstr_num = 0;
                end

            end

          end
      
      end

      updated_instance = unique( element_instance(deleted_element) );
      if length(updated_instance) > 0  && contact_flag == 1
          %'updated_instance'
          for i = updated_instance
              %i
              
              %--- Contact Pair base 
              for j = 1 : length( MODEL.CP )

                  u_ele = [];
                  for k = 1 : MODEL.INSTANCE(i).nElement 
                     offset = MODEL.INSTANCE(i).element_offset; 
                     if element_flag(k+offset) == 1
                        u_ele = [u_ele k];
                     end
                  end

                  if MODEL.CP(j).instance_id_i == i          
                      [c_triangles, c_triangles_eleid, c_nodes] = get_surface_triangle(MODEL.INSTANCE(i), u_ele, MODEL.CP(j).elements_i);
                      MODEL.CP(j).c_triangles_i = c_triangles;
                      MODEL.CP(j).c_triangles_eleid_i = c_triangles_eleid;
                      MODEL.CP(j).c_nodes_i = c_nodes;
                      %'1'
                      %c_triangles'
                      %c_triangles_eleid'
                      %c_nodes
                  end

                  if MODEL.CP(j).instance_id_j == i            
                      [c_triangles, c_triangles_eleid, c_nodes] = get_surface_triangle(MODEL.INSTANCE(i), u_ele, MODEL.CP(j).elements_j);
                      MODEL.CP(j).c_triangles_j = c_triangles;
                      MODEL.CP(j).c_triangles_eleid_j = c_triangles_eleid;
                      MODEL.CP(j).c_nodes_j = c_nodes;
                      %'2'
                      %c_triangles'
                      %c_triangles_eleid'
                      %c_nodes
                  end

              end
    
          end
      end
      

      if rem(t,d_out) == 0
        output_data(i_out).disp = disp;
        output_data(i_out).integ_stress = integ_stress;
        output_data(i_out).integ_strain = integ_strain;
        output_data(i_out).integ_plastic_strain = integ_plastic_strain;
        output_data(i_out).integ_eq_plastic_strain = integ_eq_plastic_strain;
        output_data(i_out).integ_triax_stress = integ_triax_stress;
        output_data(i_out).element_flag = element_flag;

        [node_value] = cal_node_stress_strain(nNode, elementmat, integ_num, output_data(i_out) );
        write_vtk(i_out, coordmat, elementmat, output_data(i_out).element_flag, output_data(i_out).disp, velo, node_value);

        i_out = i_out + 1;
      end

    end

    fprintf('elapsed time[s]\n')
    toc

    

    drawElement(position, elementmat, 0);

    figure;
    plot(output_);

    % for i = 1 : output_num
    %     [node_value] = cal_node_stress_strain(nNode, elementmat, integ_num, output_data(i) );
    %     write_vtk(i, coordmat, elementmat, output_data(i).element_flag, output_data(i).disp, node_value);
    % end

    res.diag_M = diag_M;
    res.diag_C = diag_C;
    res.element = elementmat;
    res.element_flag = element_flag;
    res.coordmat = coordmat;
    res.position = position;
    res.disp = disp;
    res.disp3 = reshape(res.disp, 3, nNode).';
    res.disp_pre = disp_pre;
    res.integ_stress = integ_stress;
    res.integ_strain = integ_strain;
    res.integ_yield_stress = integ_yield_stress;
    res.integ_eq_plastic_strain = integ_eq_plastic_strain;
    res.integ_triax_stress = integ_triax_stress;
    res.c_force3 = c_force3;
    res.d_node = d_node;
    
    output_data(i_out).disp = disp;
    output_data(i_out).integ_stress = integ_stress;
    output_data(i_out).integ_strain = integ_strain;
    output_data(i_out).integ_plastic_strain = integ_plastic_strain;
    output_data(i_out).integ_eq_plastic_strain = integ_eq_plastic_strain;
    output_data(i_out).integ_triax_stress = integ_triax_stress;
    output_data(i_out).element_flag = element_flag;   
    [node_value] = cal_node_stress_strain(nNode, elementmat, integ_num, output_data(i_out) );
    res.node_value = node_value;
    res.Q = Q;
    res.Pusai_mat = Pusai_mat;
    res.MODEL = MODEL;

    fclose(bug_report);
  end


  function t = cal_triax_stress_e( s )

        t = 0;
        ox = s(1);
        oy = s(2);
        oz = s(3);
        txy = s(4);
        tyz = s(5);
        txz = s(6);

        T = [ox  txy txz
             txy oy  tyz
             txz tyz oz ];
        p = eig(T);
        %v = (p(1)+p(2)+p(3))/3 / oeq;

        oeq = sqrt(0.5*( (p(1)-p(2))^2 + (p(2)-p(3))^2 + (p(3)-p(1))^2  ));
        if oeq < 1E-10  %oeq == 0
            return
        end

        t = (p(1)+p(2)+p(3))/3 / oeq;
        
  end


  function integ_triax_stress = cal_triax_stress( integ_stress )
    %n = size(integ_stress,1);
    n = size(integ_stress,2);
    integ_triax_stress = zeros(n,1);

    for i =  1 : n
        ox = integ_stress(1,i);
        oy = integ_stress(2,i);
        oz = integ_stress(3,i);
        txy = integ_stress(4,i);
        tyz = integ_stress(5,i);
        txz = integ_stress(6,i);

        %oeq = sqrt(0.5*( (ox-oy)^2 + (oy-oz)^2 + (ox-oz)^2 + 6*(txy^2 + tyz^2 + txz^2) ));
        % if oeq == 0
        %     continue
        % end

        T = [ox  txy txz
             txy oy  tyz
             txz tyz oz ];
        p = eig(T);
        %v = (p(1)+p(2)+p(3))/3 / oeq;

        oeq = sqrt(0.5*( (p(1)-p(2))^2 + (p(2)-p(3))^2 + (p(3)-p(1))^2  ));
        if oeq < 1E-10  %oeq == 0
            continue
        end

        v = (p(1)+p(2)+p(3))/3 / oeq;
        integ_triax_stress(i) = v;
    end

  end
  

  function [d_integ_stress_, d_integ_strain_, d_integ_yield_stress_, d_integ_plastic_strain_, d_integ_eq_plastic_strain_, Q] =...
            cal_stress_hexa(position_, d_disp_, elementmat, element_flag, integ_num, Pusai_mat, integ_stress_pre_, ...
                          MATERIAL, element_material, integ_yield_stress_, integ_eq_plastic_strain_, elementMinSize)

    nElement = size(elementmat,2); %size(elementmat,1);
    d_integ_stress_ = zeros(6, nElement*integ_num);
    d_integ_strain_ = zeros(6, nElement*integ_num);
    d_integ_plastic_strain_ = zeros(6, nElement*integ_num);
    d_integ_yield_stress_ = zeros(nElement*integ_num, 1);
    d_integ_eq_plastic_strain_ = zeros(nElement*integ_num, 1);
    
    fn = size(position_,2)*3;
    Q = zeros(fn, 1);
    
     % weight = 1
     W = [2];
     if integ_num == 8
        W = ones(8,1);
     end

     % mat_id = element_material(1);
     % Dmat = MATERIAL(mat_id).Dmat;
     % G = MATERIAL(mat_id).G;
     % plastic_property_ = MATERIAL(mat_id).plastic; %->even

     for e = 1 : nElement

         if element_flag(e) == 0
            continue;
         end

         mat_id = element_material(e);
         Dmat = MATERIAL(mat_id).Dmat;
         G = MATERIAL(mat_id).G;
         plastic_property_ = MATERIAL(mat_id).plastic;

        e_position = zeros(3,8);

        d_u = zeros(24,1);
        for i = 1 : 8
            e_position(1,i) = position_(1, elementmat(i,e));
            e_position(2,i) = position_(2, elementmat(i,e));
            e_position(3,i) = position_(3, elementmat(i,e)); % -> even

            %d_u((1:3)+(i-1)*3) = d_disp_( (1:3)+(elementmat(e,i)-1)*3 ); %-> slow
            d_u(1+(i-1)*3) = d_disp_( 1+(elementmat(i,e)-1)*3 );
            d_u(2+(i-1)*3) = d_disp_( 2+(elementmat(i,e)-1)*3 );
            d_u(3+(i-1)*3) = d_disp_( 3+(elementmat(i,e)-1)*3 );
        end
        
        % if norm(d_u) < elementMinSize * 1.0E-7
        %     %norm(d_u)
        %     continue
        % end

        q_vec = zeros(24,1);

        Barray = zeros(6,24,integ_num);
        %Jarray(1:integ_num) = struct('mat',zeros(3, 3) ); %->slow
        %Jarray = cell(integ_num,1); %->slow
        detJarray = zeros(integ_num,1);
        
        BVarray = zeros(6,24,integ_num);
        BVbar = zeros(6,24);
        V = 0;
        for i = 1 : integ_num
            [B, J, P2] = cal_B_hexa(Pusai_mat(:,:,i), e_position);
            Barray(:,:,i) = B;
            detJi = my3det(J);
            detJarray(i) = detJi;
            V = V + detJi;
            
            [BV_i, BVbar_i] = cal_BVbar( P2, detJi );
            BVarray(:,:,i) = BV_i;
            BVbar = BVbar + BVbar_i;
        end
        BVbar = BVbar / V;
        %V
        %Barray
        %Jarray

        
        pre_stress = zeros(6,1);
        tri_dev_stress = zeros(6,1);
        mean_stress_vec = zeros(6,1);
        for i = 1 : integ_num
          %B = Barray(:,:,i);
          B = Barray(:,:,i) + BVbar - BVarray(:,:,i);

          d_e_vec = B * d_u;
          d_o_vec = Dmat * d_e_vec;

          index_i = (e-1)*integ_num+i;
          %pre_stress = integ_stress_pre_(:,index_i);
          pre_stress(1) = integ_stress_pre_(1,index_i);
          pre_stress(2) = integ_stress_pre_(2,index_i);
          pre_stress(3) = integ_stress_pre_(3,index_i);
          pre_stress(4) = integ_stress_pre_(4,index_i);
          pre_stress(5) = integ_stress_pre_(5,index_i);
          pre_stress(6) = integ_stress_pre_(6,index_i);

          if ~isempty( plastic_property_ )    %length( plastic_property_ ) > 0 
               tri_stress = pre_stress + d_o_vec;
               mean_stress = ( tri_stress(1)+tri_stress(2)+tri_stress(3) ) / 3;
               % tri_dev_stress = [tri_stress(1) - mean_stress
               %                   tri_stress(2) - mean_stress
               %                   tri_stress(3) - mean_stress
               %                   tri_stress(4)
               %                   tri_stress(5)
               %                   tri_stress(6)];
               tri_dev_stress(1) = tri_stress(1) - mean_stress;  % -> fast
               tri_dev_stress(2) = tri_stress(2) - mean_stress;
               tri_dev_stress(3) = tri_stress(3) - mean_stress;
               tri_dev_stress(4) = tri_stress(4);
               tri_dev_stress(5) = tri_stress(5);
               tri_dev_stress(6) = tri_stress(6);
               % tri_mises_stress = sqrt( 3/2 * (tri_dev_stress(1)^2 + tri_dev_stress(2)^2 + tri_dev_stress(3)^2 + ...
               %                                 2*tri_dev_stress(4)^2 + 2*tri_dev_stress(5)^2 + 2*tri_dev_stress(6)^2) );
               tri_mises_stress = sqrt( 1.5 * (tri_dev_stress(1) * tri_dev_stress(1) + ...
                                               tri_dev_stress(2) * tri_dev_stress(2) + ...
                                               tri_dev_stress(3) * tri_dev_stress(3) + ...
                                               2*tri_dev_stress(4) * tri_dev_stress(4) + ...
                                               2*tri_dev_stress(5) * tri_dev_stress(5) + ...
                                               2*tri_dev_stress(6) * tri_dev_stress(6) ) );  %->even
               y = integ_yield_stress_( index_i );
               if tri_mises_stress > y
                   p_index = 1;
                   for j = 2 : size(plastic_property_,1)
                       if integ_eq_plastic_strain_( index_i ) <= plastic_property_(j,2)
                            p_index = j-1;
                            break
                       end
                       if j == size(plastic_property_,1)
                            p_index = j-1;
                       end
                   end
                   %p_index
                   H = (plastic_property_(p_index+1,1) - plastic_property_(p_index,1)) / (plastic_property_(p_index+1,2) - plastic_property_(p_index,2));
                   d_ep = ( tri_mises_stress - y ) / (3*G + H);
                   d_integ_eq_plastic_strain_( index_i ) = d_ep;
                   d_integ_yield_stress_( index_i ) = H * d_ep;
                   final_dev_stress = tri_dev_stress * (y+H*d_ep) / tri_mises_stress;
                   %final_stress = final_dev_stress + [mean_stress mean_stress mean_stress 0 0 0].';
                   mean_stress_vec(1) = mean_stress;
                   mean_stress_vec(2) = mean_stress;
                   mean_stress_vec(3) = mean_stress;
                   final_stress = final_dev_stress + mean_stress_vec;
                   %final_stress.'
                   %d_o_vec_tri = d_o_vec 
                   %d_o_vec = final_stress - pre_stress.';
                   d_o_vec = final_stress - pre_stress;
               end

          end

          %d_integ_stress_(:, index_i) = d_o_vec;
		  %d_integ_strain_(:, index_i) = d_e_vec;
          d_integ_stress_(1, index_i) = d_o_vec(1);
          d_integ_stress_(2, index_i) = d_o_vec(2);
          d_integ_stress_(3, index_i) = d_o_vec(3);
          d_integ_stress_(4, index_i) = d_o_vec(4);
          d_integ_stress_(5, index_i) = d_o_vec(5);
          d_integ_stress_(6, index_i) = d_o_vec(6);

          d_integ_strain_(1, index_i) = d_e_vec(1);
          d_integ_strain_(2, index_i) = d_e_vec(2);
          d_integ_strain_(3, index_i) = d_e_vec(3);
          d_integ_strain_(4, index_i) = d_e_vec(4);
          d_integ_strain_(5, index_i) = d_e_vec(5);
          d_integ_strain_(6, index_i) = d_e_vec(6);

          %o_vec = pre_stress.' + d_o_vec;
          o_vec = pre_stress + d_o_vec;
          q_vec_i = B.' * o_vec;
          q_vec = q_vec + W(i)*W(i)*W(i) * detJarray(i) * q_vec_i;

        end %i = 1 : integ_num

        for i = 1 : 8
            Q( 1 + (elementmat(i,e)-1)*3 ) = Q( 1 + (elementmat(i,e)-1)*3 ) + q_vec(1 + (i-1)*3 );
            Q( 2 + (elementmat(i,e)-1)*3 ) = Q( 2 + (elementmat(i,e)-1)*3 ) + q_vec(2 + (i-1)*3 );
            Q( 3 + (elementmat(i,e)-1)*3 ) = Q( 3 + (elementmat(i,e)-1)*3 ) + q_vec(3 + (i-1)*3 );
        end
        %Q( elementdofmat(:,e) ) = Q( elementdofmat(:,e) ) + q_vec; %->slow

     end

  end


  function [BV_i, BVbar_i] = cal_BVbar(P2, detJ_i)

    BV_i = zeros(6,24);
    BVbar_i = zeros(6,24);

    N = reshape(P2,1,24);

    % BV_i(1,:) = N;
    % BV_i(2,:) = N;
    % BV_i(3,:) = N;
    for i = 1 : 24
        BV_i(1,i) = N(i);
        BV_i(2,i) = N(i);
        BV_i(3,i) = N(i);
    end
    BV_i = BV_i / 3;
    
    BVbar_i = BV_i * detJ_i;

  end


  function [B, J, P2] = cal_B_hexa(Pusai1, e_position )

    J = Pusai1 * e_position';
    P2 = my3inv(J) * Pusai1; 
    
    B = zeros(6,24);
    
    for i = 1 : 8
      B(1,(i-1)*3+1) = P2(1,i);
      B(2,(i-1)*3+2) = P2(2,i);
      B(3,(i-1)*3+3) = P2(3,i);
      B(4,(i-1)*3+1) = P2(2,i);
      B(4,(i-1)*3+2) = P2(1,i);
      B(5,(i-1)*3+2) = P2(3,i);
      B(5,(i-1)*3+3) = P2(2,i);
      B(6,(i-1)*3+1) = P2(3,i);
      B(6,(i-1)*3+3) = P2(1,i);
    end

  end


  function Pusai_mat = cal_Pusai_hexa(integ_num)

    Pusai_mat = zeros(3,8,integ_num);
     
    delta_mat = [ -1.0  -1.0  -1.0
                   1.0  -1.0  -1.0
                   1.0   1.0  -1.0
                  -1.0   1.0  -1.0
                  -1.0  -1.0   1.0
                   1.0  -1.0   1.0
                   1.0   1.0   1.0
                  -1.0   1.0   1.0 ];

     % integral point
     gc = [0 0 0]; 

     if integ_num == 8
         g = 1.0 / sqrt(3);
         gc = [-g -g -g
               -g -g  g
               -g  g -g
               -g  g  g
                g -g -g
                g -g  g
                g  g -g
                g  g  g];
     end

    for k = 1 : integ_num
     Pusai1 = zeros(3,8);
     gzai = gc(k,1);
     eta = gc(k,2);
     tueta = gc(k,3);

      for i = 1 : 8
        Pusai1(1,i) = 1.0/8.0*delta_mat(i,1)*(1.0+eta* delta_mat(i,2))*(1.0+tueta* delta_mat(i,3));
        Pusai1(2,i) = 1.0/8.0*delta_mat(i,2)*(1.0+gzai* delta_mat(i,1))*(1.0+tueta* delta_mat(i,3));
        Pusai1(3,i) = 1.0/8.0*delta_mat(i,3)*(1.0+gzai* delta_mat(i,1))*(1.0+eta* delta_mat(i,2));
      end

      Pusai_mat(:,:,k) = Pusai1;

    end


  end


  function [element_ctr] = cal_element_ctr(position, elementmat)

        nElement = size(elementmat,2); %size(elementmat,1);
        element_ctr = zeros(3, nElement); %zeros(nElement,3);

        for i = 1 : nElement
            element_ctr(:,i) = sum(position(:,elementmat(:,i)),2) / 8;
        end

  end


  function [faces, faces_eleid, sorted_faces] = get_element_face(MODEL, i)

        part_id = MODEL.INSTANCE(i).part_id;
        cdmat = MODEL.PART( part_id ).coordmat;
        nE = MODEL.INSTANCE(i).nElement;
        faces = zeros(nE*6, 4);
        faces_eleid = zeros(nE*6, 1);
        sorted_faces = zeros(nE*6, 4);

        c = 1;
        for j = 1 : nE
            elem = MODEL.PART( part_id ).elementmat(:,j);
            faces(6*(c-1)+1, :) = [elem(1:4)];
            faces(6*(c-1)+2, :) = [elem(5:8)];
            faces(6*(c-1)+3, :) = [elem(1) elem(2) elem(6) elem(5)];
            faces(6*(c-1)+4, :) = [elem(2) elem(3) elem(7) elem(6)];
            faces(6*(c-1)+5, :) = [elem(3) elem(4) elem(8) elem(7)];
            faces(6*(c-1)+6, :) = [elem(4) elem(1) elem(5) elem(8)];
            faces_eleid(6*(c-1)+(1:6)) = j;
            ctr = sum( cdmat(:,elem),2) / 8;
            for k = 1 : 6
                index = 6*(c-1)+k;
                v1 = cdmat( :, faces(index, 2)) - cdmat( :, faces(index, 1));
                v2 = cdmat( :, faces(index, 4)) - cdmat( :, faces(index, 1));
                nv = cross(v1,v2);
                vc = ctr - cdmat( :, faces(index, 1));
                if dot(nv, vc) > 0  %nv * vc.' > 0
                    %'order modified'
                    faces(index, :) = [faces(index, 1)  faces(index, 4) faces(index, 3) faces(index, 2)];
                end
            end
            c = c + 1;
        end

        for j = 1 : nE*6
            sorted_faces(j,:) = sort(faces(j, :));
        end

  end


  function [c_triangles, c_triangles_eleid, c_nodes] = get_surface_triangle(INSTANCE_i, array_element, contact_element)

        nE = length(array_element); %MODEL.PART( part_id ).nElement;
        surfaces = zeros(nE*6, 4);
        surfaces_eleid = zeros(nE*6, 1);
        sorted_surfaces = zeros(nE*6, 4);
        c = 1;
        for j = array_element
            for k = 1 : 6
                surfaces(6*(c-1)+k, :) = INSTANCE_i.surfaces(6*(j-1)+k,:);
                sorted_surfaces(6*(c-1)+k, :) = INSTANCE_i.sorted_surfaces(6*(j-1)+k,:);
            end
            surfaces_eleid(6*(c-1)+(1:6)) = j;
            c = c + 1;
        end


        c_surfaces = [];
        c_surfaces_eleid = [];
        dp_id = [];
        for j = 1 : nE*6-1
            u_flag = 1;
            if length(find(dp_id==j)) > 0
                u_flag = 0;
                continue;
            end

            sj = sorted_surfaces(j, :);
            
            for k = j+1 : nE*6  %1 : nE*6
                % if j == k
                %     continue
                % end

                sk = sorted_surfaces(k, :);
                if sj(1) == sk(1) && sj(2) == sk(2) && sj(3) == sk(3) && sj(4) == sk(4)
                    u_flag = 0;
                    dp_id = [dp_id k];
                    break
                end
            end

            if u_flag == 1
                c_surfaces = [c_surfaces; surfaces(j, :)];
                c_surfaces_eleid = [c_surfaces_eleid  surfaces_eleid(j)];
            end

            % mat1 = sorted_surfaces(j+1 : nE*6, :);
            % mat2 = sum(abs(mat1 - ones(size(mat1,1),1) * sj).');
            % v = min(mat2);  % -> slow
            % 
            % if v > 0
            %     c_surfaces = [c_surfaces; surfaces(j, :)];
            %     c_surfaces_eleid = [c_surfaces_eleid  surfaces_eleid(j)];
            % end
        end
        length(c_surfaces);
        %MODEL.INSTANCE(i).c_surfaces = c_surfaces;
        %MODEL.INSTANCE(i).c_surfaces_eleid = c_surfaces_eleid;

        %--- Pick up only contact element
        if INSTANCE_i.nElement ~= length(contact_element)
            c_surfaces_temp = [];
            c_surfaces_eleid_temp = [];
            for j = 1 : length(c_surfaces_eleid)
                if length( find(contact_element==c_surfaces_eleid(j)) ) > 0
                    c_surfaces_temp = [c_surfaces_temp; c_surfaces(j, :)];
                    c_surfaces_eleid_temp = [c_surfaces_eleid_temp  c_surfaces_eleid(j)];
                end
            end

            c_surfaces = c_surfaces_temp;
            c_surfaces_eleid = c_surfaces_eleid_temp;
        end

        if length(c_surfaces) == 0
            c_triangles = [];
            c_triangles_eleid = [];
            c_nodes = [];
            return
        end

        c_triangles = zeros(size(c_surfaces,1)*2,3);
        c_triangles_eleid = zeros(size(c_surfaces,1)*2,1);
        for j = 1 : size(c_surfaces,1)
            c_triangles(j*2-1,:) = [c_surfaces(j,1) c_surfaces(j,2) c_surfaces(j,3)];
            c_triangles(j*2,:)   = [c_surfaces(j,3) c_surfaces(j,4) c_surfaces(j,1)];
            c_triangles_eleid(j*2-1) = c_surfaces_eleid(j);
            c_triangles_eleid(j*2) = c_surfaces_eleid(j);
        end
        
        c_nodes = reshape(c_surfaces, 1, size(c_surfaces,1)*4);
        c_nodes = unique(c_nodes);
        %length(c_nodes);
        
  end


  function [c_force3, d_node] = cal_contact_force(CP, position, velo, diag_M, elementMinSize, elementMaxSize, d_max, d_node_pre, ...
                                        INSTANCE, MATERIAL, element_flag, elementmat, bug_report, time_)  %PART, 

    c_force3 = zeros(size(position));   
    d_node = zeros(length(velo)/3, 1);
    d_lim = elementMinSize * 0.3;
    myu = 0.25;
    kc = 1.0;
    kc_s = 1.0; % self-contact

    %'cal_contact_force'

    
    instance_pair = [];
    cp_index = [];

    for cc = 1 : length(CP)

        if CP(cc).instance_id_i == CP(cc).instance_id_j
            instance_pair = [instance_pair; CP(cc).instance_id_i  CP(cc).instance_id_j];
            cp_index = [cp_index cc];
        else
            instance_pair = [instance_pair; CP(cc).instance_id_i  CP(cc).instance_id_j];
            instance_pair = [instance_pair; CP(cc).instance_id_j  CP(cc).instance_id_i];
            cp_index = [cp_index cc cc];
        end

    end

    for c = 1 : length(cp_index)

        cc = cp_index(c);

        % j -> triangle,   i -> point
        i_instance = instance_pair(c,1);  %CP(cc).instance_id_i;
        j_instance = instance_pair(c,2); %CP(cc).instance_id_j;
        if CP(cc).instance_id_i == i_instance
            c_nodes_i = CP(cc).c_nodes_i; 
            c_nodes_j = CP(cc).c_nodes_j; 
            c_triangles = CP(cc).c_triangles_j; 
            c_triangles_eleid = CP(cc).c_triangles_eleid_j;
        else
            c_nodes_i = CP(cc).c_nodes_j; 
            c_nodes_j = CP(cc).c_nodes_i; 
            c_triangles = CP(cc).c_triangles_i; 
            c_triangles_eleid = CP(cc).c_triangles_eleid_i;
        end

        node_offset_i = INSTANCE(i_instance).node_offset;
        node_offset_j = INSTANCE(j_instance).node_offset;
        element_offset_j = INSTANCE(j_instance).element_offset;
        young = MATERIAL( INSTANCE(j_instance).material_id ).young;

        c_nodes_i = c_nodes_i + node_offset_i;
        c_nodes_j = c_nodes_j + node_offset_j;
        c_triangles = c_triangles + node_offset_j;
        c_triangles_eleid = c_triangles_eleid + element_offset_j;

        %c
        %c_nodes_i
        %c_nodes_j
        %c_triangles_eleid'
        %c_triangles'

        %--- contact range ---%        

        if i_instance == j_instance % && max(d_node_pre) == 0

            %tic
                %nodeNormal = zeros(3, length(c_nodes_i));
                triangleNormal = zeros(3, size(c_triangles,1));
                triangle_v1 = zeros(3,size(c_triangles,1));
                triangle_v2 = zeros(3,size(c_triangles,1));
                triangle_q0 = zeros(3,size(c_triangles,1));
                ii = zeros(3,1);
                v1 = zeros(3,1);
                v2 = zeros(3,1);
                q0 = zeros(3,1);
                %oi = node_offset_i;
                for i = 1 : size(c_triangles,1)
                    % ii = c_triangles(i,:);
                    % q0 = position(:,ii(1));
                    % v1 = position(:,ii(2)) - q0;
                    % v2 = position(:,ii(3)) - q0;

                    %triangle_q0(:,i) = q0; 
                    % --> cost, slower than scalar input
                    %triangle_q0(i,:) = q0'; % --> slower than column major

                    % for j = 1 : 3   % faster than vector input
                    %     triangle_q0(j,i) = q0(j);
                    % end

                    ii(1) = c_triangles(i,1);
                    ii(2) = c_triangles(i,2);
                    ii(3) = c_triangles(i,3);
                    q0(1) = position(1,ii(1));
                    q0(2) = position(2,ii(1));
                    q0(3) = position(3,ii(1));
                    v1(1) = position(1,ii(2)) - q0(1);
                    v1(2) = position(2,ii(2)) - q0(2);
                    v1(3) = position(3,ii(2)) - q0(3);
                    v2(1) = position(1,ii(3)) - q0(1);
                    v2(2) = position(2,ii(3)) - q0(2);
                    v2(3) = position(3,ii(3)) - q0(3);

                    triangle_q0(1,i) = q0(1);
                    triangle_q0(2,i) = q0(2);
                    triangle_q0(3,i) = q0(3);
                    triangle_v1(1,i) = v1(1);
                    triangle_v1(2,i) = v1(2);
                    triangle_v1(3,i) = v1(3);
                    triangle_v2(1,i) = v2(1);
                    triangle_v2(2,i) = v2(2);
                    triangle_v2(3,i) = v2(3);

                    n = my3cross(v1,v2);
                    n = n / my3norm(n);
                    triangleNormal(1,i) = n(1);
                    triangleNormal(2,i) = n(2);
                    triangleNormal(3,i) = n(3);

                    % nodeNormal(1,ii(1)-oi) = nodeNormal(1,ii(1)-oi) + n(1);
                    % nodeNormal(2,ii(1)-oi) = nodeNormal(2,ii(1)-oi) + n(2);
                    % nodeNormal(3,ii(1)-oi) = nodeNormal(3,ii(1)-oi) + n(3);
                    % nodeNormal(1,ii(2)-oi) = nodeNormal(1,ii(2)-oi) + n(1);
                    % nodeNormal(2,ii(2)-oi) = nodeNormal(2,ii(2)-oi) + n(2);
                    % nodeNormal(3,ii(2)-oi) = nodeNormal(3,ii(2)-oi) + n(3);
                    % nodeNormal(1,ii(3)-oi) = nodeNormal(1,ii(3)-oi) + n(1);
                    % nodeNormal(2,ii(3)-oi) = nodeNormal(2,ii(3)-oi) + n(2);
                    % nodeNormal(3,ii(3)-oi) = nodeNormal(3,ii(3)-oi) + n(3);

                end

                A = zeros(3,3);
                p = zeros(3,1);
                b = zeros(3,1);
                n = zeros(3,1);
                for j = 1 : size(c_triangles,1)
                    j0 = c_triangles(j,1);
                    j1 = c_triangles(j,2);
                    j2 = c_triangles(j,3);
                    q0(1) = triangle_q0(1,j);
                    q0(2) = triangle_q0(2,j);
                    q0(3) = triangle_q0(3,j);
                    v1(1) = triangle_v1(1,j);
                    v1(2) = triangle_v1(2,j);
                    v1(3) = triangle_v1(3,j);
                    v2(1) = triangle_v2(1,j);
                    v2(2) = triangle_v2(2,j);
                    v2(3) = triangle_v2(3,j);
                    n(1) = triangleNormal(1,j);
                    n(2) = triangleNormal(2,j);
                    n(3) = triangleNormal(3,j);
                    %A = [v1 v2 -n]; % huge cost
                    A(1,1) = v1(1);
                    A(2,1) = v1(2);
                    A(3,1) = v1(3);
                    A(1,2) = v2(1);
                    A(2,2) = v2(2);
                    A(3,2) = v2(3);
                    A(1,3) = -n(1);
                    A(2,3) = -n(2);
                    A(3,3) = -n(3);

                    L1 = my3norm(v1);
                    L2 = my3norm(v2);
                    Lmax = max(L1,L2);
                    d12 = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3);
                    S = 0.5 * sqrt( L1*L1 * L2*L2 - d12*d12 );

                    for i = c_nodes_i
                        if j0 == i || j1 == i || j2 == i
                            continue
                        end

                        p(1) = position(1,i);
                        p(2) = position(2,i);
                        p(3) = position(3,i);
                        b(1) = p(1) -q0(1);
                        b(2) = p(2) -q0(2);
                        b(3) = p(3) -q0(3);

                        %d = my3norm(b);
                        nd = b(1)*n(1) + b(2)*n(2) + b(3)*n(3);
                        if nd > 0  ||  nd < -d_lim %elementMaxSize
                            continue;
                        end

                        if my3norm(b) > elementMaxSize
                            continue;
                        end

                        %x = A \ b;
                        x = my3inv(A) * b;
                        d = x(3);

                        if 0.0 <= x(1) && 0.0 <= x(2) && x(1) + x(2) <= 1.0 && ...
                                d > 0 && d <= d_lim

                            if d - d_node_pre(i) > d_max
                                d = d_node_pre(i) + d_max;
                            end
                            [i d n'];

                            v(1) = velo(i*3-2) - velo(j0*3-2); 
                            v(2) = velo(i*3-1) - velo(j0*3-1);
                            v(3) = velo(i*3) - velo(j0*3);
                            mag_v = my3norm(v);
                            %cv = 1.0;
                            ve = zeros(3,1);
                            if mag_v > 0
                                ve = v / my3norm(v);
                                %cv = abs(ve(1)*n(1)+ve(2)*n(2)+ve(3)*n(3));
                            end
                            %cv
                            k = young * S / Lmax * kc_s; %1.0;
                            %damp = 2 * sqrt( diag_M(i) * k ) * 0.1;
                            F = k * d;
                            f = F * n; %  - damp * dot(v,n) * n;

                            %--- Friction
                                %vs = ve - dot(ve,n) * n;
                            %vs = ve - (ve(1)*n(1)+ve(2)*n(2)+ve(3)*n(3)) * n;
                            %fric = -myu * F * vs;
                            %f = f + fric;  

                            dot_ve_n = ve(1)*n(1)+ve(2)*n(2)+ve(3)*n(3); 
                            vs(1) = ve(1) - dot_ve_n * n(1); % ->declare as vector(3) explicitly
                            vs(2) = ve(2) - dot_ve_n * n(2);
                            vs(3) = ve(3) - dot_ve_n * n(3);
                            fric = -myu * F * vs;
                            f = f + fric;  

                            c_force3(1,i) = c_force3(1,i) + f(1); 
                            c_force3(2,i) = c_force3(2,i) + f(2); 
                            c_force3(3,i) = c_force3(3,i) + f(3); 
                            c_force3(1,j0) = c_force3(1,j0) - f(1)  / 3;
                            c_force3(2,j0) = c_force3(2,j0) - f(2)  / 3;
                            c_force3(3,j0) = c_force3(3,j0) - f(3)  / 3;
                            c_force3(1,j1) = c_force3(1,j1) - f(1)  / 3;
                            c_force3(2,j1) = c_force3(2,j1) - f(2)  / 3;
                            c_force3(3,j1) = c_force3(3,j1) - f(3)  / 3;
                            c_force3(1,j2) = c_force3(1,j2) - f(1)  / 3;
                            c_force3(2,j2) = c_force3(2,j2) - f(2)  / 3;
                            c_force3(3,j2) = c_force3(3,j2) - f(3)  / 3;


                            if d > d_node(i)
                                d_node(i) = d;
                            end

                        end


                    end
                end

                %'node normal'
            %toc

            continue;
        end


        %--- contact range ---%    
        %if i_instance ~= j_instance
            position_i = position(:, c_nodes_i);
            position_j = position(:, c_nodes_j);
            min_i = min(position_i,[],2);
            max_i = max(position_i,[],2);
            min_j = min(position_j,[],2);
            max_j = max(position_j,[],2);
    
            range_min_x = max(min_i(1),min_j(1));
            range_max_x = min(max_i(1),max_j(1));
            range_min_y = max(min_i(2),min_j(2));
            range_max_y = min(max_i(2),max_j(2));
            range_min_z = max(min_i(3),min_j(3));
            range_max_z = min(max_i(3),max_j(3));
    
            if range_min_x > range_max_x || range_min_y > range_max_y || range_min_z > range_max_z
                continue
            end
        %end
        
        A = zeros(3,3);
        p = zeros(3,1);
        q0 = zeros(3,1);
        q1 = zeros(3,1);
        q2 = zeros(3,1);
        elem = zeros(8,1);
        for j = 1 : size(c_triangles,1)
            if element_flag( c_triangles_eleid(j) ) == 0
                continue
            end
            j0 = c_triangles(j,1);
            j1 = c_triangles(j,2);
            j2 = c_triangles(j,3);
            
            q0(1) = position(1,j0);
            q0(2) = position(2,j0);
            q0(3) = position(3,j0);
            q1(1) = position(1,j1);
            q1(2) = position(2,j1);
            q1(3) = position(3,j1);
            q2(1) = position(1,j2);
            q2(2) = position(2,j2);
            q2(3) = position(3,j2);
            
            %if i_instance ~= j_instance
                if q0(1) < range_min_x && q1(1) < range_min_x && q2(1) < range_min_x
                    continue
                end
                if q0(2) < range_min_y && q1(2) < range_min_y && q2(2) < range_min_y
                    continue
                end
                if q0(3) < range_min_z && q1(3) < range_min_z && q2(3) < range_min_z
                    continue
                end
    
                if q0(1) > range_max_x && q1(1) > range_max_x && q2(1) > range_max_x
                    continue
                end
                if q0(2) > range_max_y && q1(2) > range_max_y && q2(2) > range_max_y
                    continue
                end
                if q0(3) > range_max_z && q1(3) > range_max_z && q2(3) > range_max_z
                    continue
                end
            %end
        

            v1 = q1 -q0;
            v2 = q2 -q0;
            L1 = my3norm(v1);
            L2 = my3norm(v2);
            Lmax = max(L1,L2);
            n = my3cross(v1,v2); %cross(v1,v2);
            n = n / my3norm(n);  %norm(n);
            d12 = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3);
            S = 0.5 * sqrt( L1*L1 * L2*L2 - d12*d12 );
            
            %--- custom
            % if abs(n(3)) > 0.9 || abs(n(1)) > 0.9
            %     continue
            % end

            %A = [v1 v2 -n]; % huge cost
            A(1,1) = v1(1);
            A(2,1) = v1(2);
            A(3,1) = v1(3);
            A(1,2) = v2(1);
            A(2,2) = v2(2);
            A(3,2) = v2(3);
            A(1,3) = -n(1);
            A(2,3) = -n(2);
            A(3,3) = -n(3);

            for i = c_nodes_i
                p(1) = position(1,i);
                p(2) = position(2,i);
                p(3) = position(3,i);

                %if i_instance ~= j_instance
                    if p(1) < range_min_x || p(2) < range_min_y || p(3) < range_min_z
                        continue
                    end
                    if p(1) > range_max_x || p(2) > range_max_y || p(3) > range_max_z
                        continue
                    end
                %end

                b = (p - q0);
                if my3norm(b) > Lmax
                    continue
                end


                % dA = abs(det(A));
                % if dA < 1E-10
                %     dA
                %     continue;
                % end

                % rA = abs(rcond(A));
                % if rA < 1E-3 %|| isnan(rA)
                %     rA
                %     continue;
                % end
                
                %x = A \ b;
                x = my3inv(A) * b;
                d = x(3);
                
                if 0.0 <= x(1) && 0.0 <= x(2) && x(1) + x(2) <= 1.0 && ...
                        d > 0 && d <= d_lim
                     % x
                     % i
                     % p
                     % q0
                     % q1
                     % q2
                     % n
                     % dA
                     % d
                     %[i d n']
                    
                    %d_max = 6.6E-04;
                     if d - d_node_pre(i) > d_max
                         d = d_node_pre(i) + d_max;
                     end
                    %[i d n']

                    v(1) = velo(i*3-2) - velo(j0*3-2); 
                    v(2) = velo(i*3-1) - velo(j0*3-1);
                    v(3) = velo(i*3) - velo(j0*3);
                    mag_v = my3norm(v);
                        %cv = 1.0;
                    ve = zeros(3,1);
                    vs = zeros(3,1);
                    if mag_v > 0.0   % 1.0E-10  % 0          
                            %ve = v / norm(v);
                        ve = v / mag_v;
                            %cv = abs( ve(1)*n(1)+ve(2)*n(2)+ve(3)*n(3) ); %dot(ve,n)
                    end
                        %fprintf("ve:[%.6e,%.6e,%.6e], cv=%.6e\n", ve(1),ve(2),ve(3),cv)  % ve, c -> same as julia
                    k = young * S / Lmax * kc;   %10.0;
                        %damp = 2 * sqrt( diag_M(i) * k ) * 0.1;
                    F = k * d;
                    f = F * n; %  - damp * dot(v,n) * n;

                        %fprintf(bug_report, "time=%.3e, i=%d, cv=%.6e, vi=%.6e, %.6e, %.6e,\n",...
                        %                     time_, i, cv, v(1), v(2), v(3));
                        %   -> no diff. from julia
                        %fprintf(bug_report, "time=%.3e, i=%d, vi=%.16e, %.16e, %.16e\n",...
                        %                     time_, i, v(1), v(2), v(3));
                        %   -> diff. from julia.  This might be double rounding
                        %   error. not coding bugs
                        %fprintf(bug_report, "time=%.3e, i=%d, pi=%.16e, %.16e, %.16e\n",...
                        %                     time_, i, p(1), p(2), p(3));

                    %--- Friction
                    %dot_ve_n = ve(1)*n(1)+ve(2)*n(2)+ve(3)*n(3)
                    %vs = ve - (dot_ve_n * n)
                    %vs = ve - (ve(1)*n(1)+ve(2)*n(2)+ve(3)*n(3)) * n;  % -> dimension error. "vs" becomes matrix not vector

                    dot_ve_n = ve(1)*n(1)+ve(2)*n(2)+ve(3)*n(3); 
                    vs(1) = ve(1) - dot_ve_n * n(1); % ->declare as vector(3) explicitly
                    vs(2) = ve(2) - dot_ve_n * n(2);
                    vs(3) = ve(3) - dot_ve_n * n(3);
                    fric = -myu * F * vs;
                    f = f + fric;  

                    
                    c_force3(1,i) = c_force3(1,i) + f(1); 
                    c_force3(2,i) = c_force3(2,i) + f(2); 
                    c_force3(3,i) = c_force3(3,i) + f(3); 

                    c_force3(1,j0) = c_force3(1,j0) - f(1)  / 3;
                    c_force3(2,j0) = c_force3(2,j0) - f(2)  / 3;
                    c_force3(3,j0) = c_force3(3,j0) - f(3)  / 3;

                    c_force3(1,j1) = c_force3(1,j1) - f(1)  / 3;
                    c_force3(2,j1) = c_force3(2,j1) - f(2)  / 3;
                    c_force3(3,j1) = c_force3(3,j1) - f(3)  / 3;

                    c_force3(1,j2) = c_force3(1,j2) - f(1)  / 3;
                    c_force3(2,j2) = c_force3(2,j2) - f(2)  / 3;
                    c_force3(3,j2) = c_force3(3,j2) - f(3)  / 3;
                    

                    if d > d_node(i)
                        d_node(i) = d;
                    end
                   
                    
                end

            end
        end

    end %c

  end


  function v = my3norm(v1)
    v = sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3));
  end

  function n = my3cross(a,b)
    n = zeros(3,1);
    n(1) = a(2)*b(3) - a(3)*b(2);
    n(2) = a(3)*b(1) - a(1)*b(3);
    n(3) = a(1)*b(2) - a(2)*b(1);
  end

  function v = my3det(m3)
    v =   m3(1,1)*m3(2,2)*m3(3,3) ...
        + m3(1,2)*m3(2,3)*m3(3,1) ...
        + m3(1,3)*m3(2,1)*m3(3,2) ...
        - m3(1,1)*m3(2,3)*m3(3,2) ...
        - m3(1,2)*m3(2,1)*m3(3,3) ...
        - m3(1,3)*m3(2,2)*m3(3,1);
  end

  function im = my3inv(m3)
    v =   m3(1,1)*m3(2,2)*m3(3,3) ...
        + m3(1,2)*m3(2,3)*m3(3,1) ...
        + m3(1,3)*m3(2,1)*m3(3,2) ...
        - m3(1,1)*m3(2,3)*m3(3,2) ...
        - m3(1,2)*m3(2,1)*m3(3,3) ...
        - m3(1,3)*m3(2,2)*m3(3,1);

    im = zeros(3,3);
    im(1,1) = ( m3(2,2)*m3(3,3) - m3(2,3)*m3(3,2) ) / v;
    im(2,1) = ( m3(2,3)*m3(3,1) - m3(2,1)*m3(3,3) ) / v;
    im(3,1) = ( m3(2,1)*m3(3,2) - m3(2,2)*m3(3,1) ) / v;

    im(1,2) = ( m3(1,3)*m3(3,2) - m3(1,2)*m3(3,3) ) / v;
    im(2,2) = ( m3(1,1)*m3(3,3) - m3(1,3)*m3(3,1) ) / v;
    im(3,2) = ( m3(1,2)*m3(3,1) - m3(1,1)*m3(3,2) ) / v;

    im(1,3) = ( m3(1,2)*m3(2,3) - m3(1,3)*m3(2,2) ) / v;
    im(2,3) = ( m3(1,3)*m3(2,1) - m3(1,1)*m3(2,3) ) / v;
    im(3,3) = ( m3(1,1)*m3(2,2) - m3(1,2)*m3(2,1) ) / v;

  end



  function [node_value] = cal_node_stress_strain(nNode, elementmat, integ_num, data)

    integ_stress = data.integ_stress';
    integ_strain = data.integ_strain';
    integ_eq_plastic_strain = data.integ_eq_plastic_strain;
    integ_triax_stress = data.integ_triax_stress;

    node_stress = zeros(nNode, 6);
    node_strain = zeros(nNode, 6);
    node_eq_plastic_strain = zeros(nNode,1);
    node_triax_stress = zeros(nNode,1);

    nElement = size(elementmat,2); %size(elementmat,1);
    element_stress = zeros(nElement, 6);
    element_strain = zeros(nElement, 6);
    element_eq_plastic_strain = zeros(nElement, 1);
    element_triax_stress = zeros(nElement, 1);

    for i = 1 : nElement
        if integ_num == 1
            element_stress(i,:) = integ_stress(i, :);
            element_strain(i,:) = integ_strain(i, :);
            element_eq_plastic_strain(i) = integ_eq_plastic_strain(i);
            element_triax_stress(i) = integ_triax_stress(i);
        else
            element_stress(i,:) = sum( integ_stress((1:integ_num)+(i-1)*integ_num, :)  ) / integ_num;
            element_strain(i,:) = sum( integ_strain((1:integ_num)+(i-1)*integ_num, :)  ) / integ_num;
            element_eq_plastic_strain(i) = sum( integ_eq_plastic_strain((1:integ_num)+(i-1)*integ_num)  ) / integ_num;
            element_triax_stress(i) = sum( integ_triax_stress((1:integ_num)+(i-1)*integ_num)  ) / integ_num;
        end
    end

    for i = 1 : nElement
        for k = 1 : 8
            node_stress(elementmat(k,i),:) = node_stress(elementmat(k,i),:) + element_stress(i, :);
            node_strain(elementmat(k,i),:) = node_strain(elementmat(k,i),:) + element_strain(i, :);
            node_eq_plastic_strain(elementmat(k,i)) = node_eq_plastic_strain(elementmat(k,i)) + element_eq_plastic_strain(i);
            node_triax_stress(elementmat(k,i)) = node_triax_stress(elementmat(k,i)) + element_triax_stress(i);
        end
    end

    inc_num = zeros(nNode, 1);
    for i = 1 : nElement
       inc_num(elementmat(:,i)) = inc_num(elementmat(:,i)) + 1;
    end
    %inc_num

    for i = 1 : nNode
        node_stress(i,:) = node_stress(i,:) / inc_num(i);
        node_strain(i,:) = node_strain(i,:) / inc_num(i);
        node_eq_plastic_strain(i) = node_eq_plastic_strain(i) / inc_num(i);
        node_triax_stress(i) = node_triax_stress(i) / inc_num(i);
    end

    node_mises_stress = zeros(nNode,1);
    for i = 1 : nNode
        ox = node_stress(i,1);
        oy = node_stress(i,2);
        oz = node_stress(i,3);
        txy = node_stress(i,4);
        tyz = node_stress(i,5);
        txz = node_stress(i,6);
        node_mises_stress(i) = sqrt( 0.5*( (ox-oy)^2 + (oy-oz)^2 + (ox-oz)^2 + 6 * (txy^2 + tyz^2 + txz^2) ) );
    end

    node_value.node_stress = node_stress;
    node_value.node_strain = node_strain;
    node_value.node_mises_stress = node_mises_stress;
    node_value.node_eq_plastic_strain = node_eq_plastic_strain;
    node_value.node_triax_stress = node_triax_stress;

end


function write_vtk(index_, coordmat, elementmat, element_flag, disp, velo, data)

    node_stress = data.node_stress;
    node_strain = data.node_strain;
    node_mises_stress = data.node_mises_stress;
    node_eq_plastic_strain = data.node_eq_plastic_strain;
    node_triax_stress = data.node_triax_stress;

    nNode = size(coordmat,2); %size(coordmat,1);
    nElement = size(elementmat,2); %size(elementmat,1);
    disp3 = reshape(disp,3,nNode)';
    velo3 = reshape(velo,3,nNode)';

    for i = 1 : nNode
        for j = 1 : 3
            if abs(disp3(i,j)) < 1E-16
                disp3(i,j) = 0;
            end
            if abs(velo3(i,j)) < 1E-16
                velo3(i,j) = 0;
            end
        end

        for j = 1 : 6
            if abs(node_stress(i,j)) < 1E-16
                node_stress(i,j) = 0;
            end
            if abs(node_strain(i,j)) < 1E-16
                node_strain(i,j) = 0;
            end
        end

        if abs(node_mises_stress(i)) < 1E-16
            node_mises_stress(i) = 0;
        end
        if abs(node_eq_plastic_strain(i)) < 1E-16
            node_eq_plastic_strain(i) = 0;
        end
        if abs(node_triax_stress(i)) < 1E-16
            node_triax_stress(i) = 0;
        end
    end

    if not(exist('temp','dir'))
        mkdir('temp')
    end

    fname = sprintf("temp\\file%03d.vtk", index_);
    out = fopen(fname,"w");
    fprintf(out,"# vtk DataFile Version 2.0\n");
    fprintf(out,"Test\n");
    fprintf(out,"ASCII\n");
    fprintf(out,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(out,"POINTS %d float\n", nNode);
    
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e %1.6e %1.6e\n", coordmat(1,i), coordmat(2,i), coordmat(3,i)) );
    end

    draw_element_num = sum(element_flag);

    %fprintf(out, "CELLS %d %d\n", nElement, nElement*(8+1));
    fprintf(out, "CELLS %d %d\n", draw_element_num, draw_element_num*(8+1));
    e = elementmat;
    e = e - 1;
    for i = 1 : nElement
        if element_flag(i) == 1
            fprintf(out, sprintf("8 %d %d %d %d %d %d %d %d\n", ...
                e(1,i), e(2,i), e(3,i), e(4,i), e(5,i), e(6,i), e(7,i), e(8,i)) );
        end
    end

    %fprintf(out,"CELL_TYPES %d\n", nElement);
    fprintf(out,"CELL_TYPES %d\n", draw_element_num);
    for i = 1 : nElement
        if element_flag(i) == 1
            fprintf(out, "12\n");
        end
    end

    
    fprintf(out,"POINT_DATA %d\n", nNode);
    fprintf(out,"VECTORS DISPLACEMENT float\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e %1.6e %1.6e\n", disp3(i,1), disp3(i,2), disp3(i,3)) );
    end


    % fprintf(out,"POINT_DATA %d\n", nNode);
    % fprintf(out,"VECTORS VELOCITY float\n");
    % for i = 1 : nNode
    %     fprintf(out, sprintf("%1.6e %1.6e %1.6e\n", velo3(i,1), velo3(i,2), velo3(i,3)) );
    % end   % -> cannot read by paraview


    fprintf(out,"SCALARS Vx float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", velo3(i,1) ) );
    end
    
    fprintf(out,"SCALARS Vy float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", velo3(i,2) ) );
    end

    fprintf(out,"SCALARS Vz float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", velo3(i,3) ) );
    end
    
    
    fprintf(out,"SCALARS E11 float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_strain(i,1) ) );
    end

    fprintf(out,"SCALARS E22 float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_strain(i,2) ) );
    end

    fprintf(out,"SCALARS E33 float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_strain(i,3) ) );
    end

    fprintf(out,"SCALARS E12 float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_strain(i,4) ) );
    end

    fprintf(out,"SCALARS E23 float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_strain(i,5) ) );
    end

    fprintf(out,"SCALARS E13 float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_strain(i,6) ) );
    end

    fprintf(out,"SCALARS EQ_PSTRAIN float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_eq_plastic_strain(i) ) );
    end


    fprintf(out,"SCALARS S11 float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_stress(i,1) ) );
    end

    fprintf(out,"SCALARS S22 float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_stress(i,2) ) );
    end

    fprintf(out,"SCALARS S33 float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_stress(i,3) ) );
    end

    fprintf(out,"SCALARS S12 float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_stress(i,4) ) );
    end

    fprintf(out,"SCALARS S23 float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_stress(i,5) ) );
    end

    fprintf(out,"SCALARS S13 float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_stress(i,6) ) );
    end

    fprintf(out,"SCALARS MISES_STRESS float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_mises_stress(i) ) );
    end

    fprintf(out,"SCALARS TRIAX_STRESS float 1\n");
    fprintf(out,"LOOKUP_TABLE default\n");
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e\n", node_triax_stress(i) ) );
    end
    
    fclose(out);

end

