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



    varargin
       
    [ MODEL] = readInpFile(varargin{1});
    nNode = MODEL.nNode;
    coordmat = MODEL.coordmat;
    nElement = MODEL.nElement;
    elementmat = MODEL.elementmat;
    element_material = MODEL.element_material;
    element_instance = MODEL.element_instance;
    contact_flag = MODEL.contact_flag;
    
    flag_fracture = 0;
    
    mass_scaling = MODEL.mass_scaling;
    d_time = MODEL.d_time * sqrt(mass_scaling);
    end_time = MODEL.end_time;
    time_num = end_time / d_time;
    C = 0.0E3;


    force_node = [];
    force_dof = [];
    force_v = [];
        
    c_pair = [];
    if contact_flag == 1
        ni = length(MODEL.INSTANCE);
        %ni = 2;
        for i = 1 : ni
            for j = 1 : ni
                if i ~= j
                    c_pair = [c_pair; i j];
                end
            end
        end    
        c_pair
    end
    
  
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
    
    elementVolume = zeros(nElement,1);
    for e = 1 : nElement
        e_position = coordmat(elementmat(e,:), :);
        V = 0;
        for i = 1 : integ_num
            %[B, J] = cal_B_hexa(Pusai_mat(:,:,i), e_position);
            J = Pusai_mat(:,:,i) * e_position;
            V = V + det(J);
        end
        elementVolume(e) = V;
    end
    elementVolume.';

    diag_M = zeros(fn, 1);
    for i = 1 : nElement
        density = MODEL.MATERIAL( element_material(i) ).density;
        node_mass = density * elementVolume(i) / 8;
        diag_M( (elementmat(i,:)-1)*3 + 1 ) = diag_M( (elementmat(i,:)-1)*3 + 1 ) + node_mass;
        diag_M( (elementmat(i,:)-1)*3 + 2 ) = diag_M( (elementmat(i,:)-1)*3 + 2 ) + node_mass;
        diag_M( (elementmat(i,:)-1)*3 + 3 ) = diag_M( (elementmat(i,:)-1)*3 + 3 ) + node_mass;
    end

    diag_M = diag_M * mass_scaling;
    diag_C = diag_M * C;


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
    if contact_flag == 1
        tic;
    
        for i = 1 : length(MODEL.INSTANCE)  
    
            [faces, faces_eleid, sorted_faces] = get_element_face(MODEL, i);
            MODEL.INSTANCE(i).surfaces = faces;
            MODEL.INSTANCE(i).surfaces_eleid = faces_eleid;
            MODEL.INSTANCE(i).sorted_surfaces = sorted_faces;
    
            array_element = [1 : MODEL.INSTANCE(i).nElement];
            [c_triangles, c_triangles_eleid, c_nodes] = get_surface_triangle(MODEL, i, array_element);
            length(c_triangles)
            MODEL.INSTANCE(i).c_triangles = c_triangles;
            MODEL.INSTANCE(i).c_triangles_eleid = c_triangles_eleid;
            MODEL.INSTANCE(i).c_nodes = c_nodes;
    
        end
    
        'Contact set'
        toc
    end


    elementSize = zeros(nElement,3);
    for e = 1 : nElement
        e_pos = coordmat(elementmat(e,:), :);
        L1 = norm( e_pos(1,:) - e_pos(2,:) );
        L2 = norm( e_pos(1,:) - e_pos(4,:) );
        L3 = norm( e_pos(1,:) - e_pos(5,:) );
        elementSize(e,:) = [L1 L2 L3];
    end
    elementMinSize = min(min(elementSize))


    %--- Variable ---%
    Q = zeros(fn, 1);
    integ_stress = zeros( nElement * integ_num, 6);
    integ_strain = zeros( nElement * integ_num, 6);
    integ_plastic_strain = zeros( nElement * integ_num, 6);
    integ_eq_plastic_strain = zeros( nElement * integ_num, 1);
    integ_triax_stress = zeros( nElement * integ_num, 1);
    element_flag = ones(nElement, 1);
    integ_flag = ones(nElement, integ_num);
  
    integ_yield_stress = zeros( nElement * integ_num, 1);
    for i = 1 : nElement
        pp = MODEL.MATERIAL(element_material(i)).plastic;
        if length( pp ) >0
            integ_yield_stress([1:integ_num]+(i-1)*integ_num) = pp(1,1);
        end
    end

    element_ctr = zeros(nElement,3);
    element_ctr = cal_element_ctr(coordmat, elementmat);


    output_num = 40
    d_out = floor(time_num / output_num)
    output_data(1:output_num) = struct('disp',zeros(fn, 1),...
                               'integ_stress',zeros( nElement * integ_num, 6),...
                               'integ_strain',zeros( nElement * integ_num, 6),...
                               'integ_plastic_strain', zeros( nElement * integ_num, 6),...
                               'integ_eq_plastic_strain', zeros( nElement * integ_num, 1),...
                               'integ_triax_stress', zeros( nElement * integ_num, 1),...
                               'element_flag',ones(nElement, 1) );

    [node_value] = cal_node_stress_strain(nNode, elementmat, integ_num, output_data(1));
    write_vtk(0, coordmat, elementmat, output_data(1).element_flag, output_data(1).disp, node_value);

    i_out = 1;
    clstr_num = 0;

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
      
      if contact_flag == 1
          c_force3 = cal_contact_force(c_pair, position, velo, diag_M, elementMinSize, ...
              MODEL.INSTANCE, MODEL.PART, MODEL.MATERIAL, element_flag);
          external_force = external_force + reshape(c_force3.', fn,1);
      end


        % update position
        disp_new = 1.0 ./ (diag_M/d_time^2.0 + diag_C/2.0/d_time) .* ( external_force - Q + diag_M/d_time^2.0 .* (2*disp - disp_pre) + diag_C/2.0/d_time.*disp_pre );
        
        % Boundary Conditions
        %disp_new( fix_dof ) = 0;
        %disp_new( input_disp_dof ) = input_disp_mag * t / time_num ;
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
                disp_new( dof ) = v * amp;
            end
        end

        d_disp = disp_new - disp;
        disp_pre = disp;
        disp = disp_new;
        velo = d_disp / d_time;

        position = coordmat + reshape(disp,3,nNode).';
        element_ctr = cal_element_ctr(position, elementmat);


      %cal element stress, strain
      [d_integ_stress, d_integ_strain, d_integ_yield_stress, d_integ_plastic_strain, d_integ_eq_plastic_strain, Q] = ...
          cal_stress_hexa(position, d_disp, elementmat, element_flag, integ_num, Pusai_mat, integ_stress, ...
                          MODEL.MATERIAL, element_material, integ_yield_stress, integ_eq_plastic_strain, elementMinSize);

      
      integ_stress = integ_stress + d_integ_stress;
      integ_strain = integ_strain + d_integ_strain;
      integ_plastic_strain = integ_plastic_strain + d_integ_plastic_strain;
      integ_eq_plastic_strain = integ_eq_plastic_strain + d_integ_eq_plastic_strain;
      integ_yield_stress = integ_yield_stress + d_integ_yield_stress;
      integ_triax_stress = cal_triax_stress( integ_stress );

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
            %fracture_stress = MODEL.MATERIAL(mat_id).failure;
            
            nd =size(MODEL.MATERIAL(mat_id).ductile, 1);
            if nd > 0

                v_e = sum( integ_eq_plastic_strain([1:integ_num]+(i-1)*integ_num) ) / integ_num;
                t_e = sum( integ_triax_stress([1:integ_num]+(i-1)*integ_num) ) / integ_num;
                if t_e < 0
                    continue
                end

                % fr_e = 0.5;
                % if t_e <= 0.4
                %     fr_e = 1.0 + (0.5 - 1.0) / (0.4 - 0.0) * t_e;
                % end

                ductile_ = MODEL.MATERIAL(mat_id).ductile;
                fr_e = ductile_(nd,1);
                for j = 1 : nd-1
                    if t_e >= ductile_(j,2) && t_e < ductile_(j+1, 2)
                        fr_e = ductile_(j,1) + (ductile_(j+1,1) - ductile_(j,1)) / (ductile_(j+1,2) - ductile_(j,2)) * (t_e - ductile_(j,2));
                        break
                    end
                end


                %if sum(integ_flag(i,:)) == 0  &&  element_flag(i) == 1
                if v_e >= fr_e  &&  element_flag(i) == 1
                    % i
                    % v_e
                    % t_e
                    % fr_e
                    element_flag( i ) = 0;
                    integ_stress( [1:integ_num]+(i-1)*integ_num, :) = 0;
                    integ_strain( [1:integ_num]+(i-1)*integ_num, :) = 0;
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
              % d_ele = [];
              % for j = deleted_element
              %   if element_instance(j) == i
              %       d_ele = [d_ele j];
              %   end
              % end
              % 
              % u_ele = [];
              % for j = d_ele
              %   L = max(elementSize(j,:)) * 2;
              %   for k = 1 : MODEL.INSTANCE(i).nElement
              %       e = k + MODEL.INSTANCE(i).element_offset;
              %       if  j == e
              %           continue
              %       end
              %       d = norm( element_ctr(e,:) - element_ctr(j,:) );
              %       if d <= L
              %           u_ele = [u_ele k];
              %       end
              %   end
              % end
              % u_ele

              u_ele = [];
              for k = 1 : MODEL.INSTANCE(i).nElement
                 offset = MODEL.INSTANCE(i).element_offset;
                 if element_flag(k+offset) == 1
                    u_ele = [u_ele k];
                 end
              end
              %[length(u_ele) MODEL.INSTANCE(i).nElement]
    
            [c_triangles, c_triangles_eleid, c_nodes] = get_surface_triangle(MODEL, i, u_ele);
            %c_triangles
            %length(c_triangles)
            MODEL.INSTANCE(i).c_triangles = c_triangles;
            MODEL.INSTANCE(i).c_triangles_eleid = c_triangles_eleid;
            MODEL.INSTANCE(i).c_nodes = c_nodes;
    
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
        i_out = i_out + 1;
      end
    end

    'elapsed time[s]'
    toc

    drawElement(position, elementmat, 0);

    for i = 1 : output_num
        [node_value] = cal_node_stress_strain(nNode, elementmat, integ_num, output_data(i) );
        write_vtk(i, coordmat, elementmat, output_data(i).element_flag, output_data(i).disp, node_value);
    end

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
    
    % res.node_stress = node_stress;
    % res.node_strain = node_strain;
    % res.node_mises_stress = node_mises_stress;
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

  end


  function integ_triax_stress = cal_triax_stress( integ_stress )
    n = size(integ_stress,1);
    integ_triax_stress = zeros(n,1);

    for i =  1 : n
        ox = integ_stress(i,1);
        oy = integ_stress(i,2);
        oz = integ_stress(i,3);
        txy = integ_stress(i,4);
        tyz = integ_stress(i,5);
        txz = integ_stress(i,6);
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

    nElement = size(elementmat,1);
    d_integ_stress_ = zeros(nElement*integ_num, 6);
    d_integ_strain_ = zeros(nElement*integ_num, 6);
    d_integ_yield_stress_ = zeros(nElement*integ_num, 1);
    d_integ_plastic_strain_ = zeros(nElement*integ_num, 6);
    d_integ_eq_plastic_strain_ = zeros(nElement*integ_num, 1);
    %Q_element_ = zeros(nElement, 24);

    fn = size(position_,1)*3;
    Q = zeros(fn, 1);
    
     % weight = 1
     W = [2];
     if integ_num == 8
        W = ones(8,1);
     end

     for e = 1 : nElement
         if element_flag(e) == 0
            continue;
         end

         mat_id = element_material(e);
         Dmat = MATERIAL(mat_id).Dmat;
         G = MATERIAL(mat_id).G;
         plastic_property_ = MATERIAL(mat_id).plastic;

        e_position = position_(elementmat(e,:), :);

        d_u = zeros(24,1);
        for i = 1 : 8
            %d_u([1:3]+(i-1)*3) = d_disp_( [1:3]+(elementmat(e,i)-1)*3 );
            d_u(1+(i-1)*3) = d_disp_( 1+(elementmat(e,i)-1)*3 );
            d_u(2+(i-1)*3) = d_disp_( 2+(elementmat(e,i)-1)*3 );
            d_u(3+(i-1)*3) = d_disp_( 3+(elementmat(e,i)-1)*3 );
        end

        % if norm(d_u) < elementMinSize * 1.0E-7
        %     %norm(d_u)
        %     continue
        % end

        q_vec = zeros(24,1);

        Barray = zeros(6,24,integ_num);
        Jarray = zeros(3,3,integ_num);
        detJarray = zeros(integ_num);
        V = 0;
        for i = 1 : integ_num
            [B, J] = cal_B_hexa(Pusai_mat(:,:,i), e_position);
            Barray(:,:,i) = B;
            Jarray(:,:,i) = J;
            detJi = det(J);
            detJarray(i) = detJi;
            %V = V + det(J);
            V = V + detJi;
        end
        %V
        %Barray
        %Jarray

        BVarray = zeros(6,24,integ_num);
        BVbar = zeros(6,24);
        for i = 1 : integ_num
            [BV_i, BVbar_i] = cal_BVbar(Pusai_mat(:,:,i), Jarray(:,:,i), detJarray(i) );
            BVarray(:,:,i) = BV_i;
            BVbar = BVbar + BVbar_i;
        end
        BVbar = BVbar / V;


        for i = 1 : integ_num
          % B(6,24)
          %[B, J] = cal_B_hexa(Pusai_mat(:,:,i), e_position);
          %B = Barray(:,:,i);
          B = Barray(:,:,i) + BVbar - BVarray(:,:,i);
          %J = Jarray(:,:,i);

          d_e_vec = B * d_u;
          d_o_vec = Dmat * d_e_vec;

          index_i = (e-1)*integ_num+i;
          pre_stress = integ_stress_pre_(index_i,:);

          if length( plastic_property_ ) > 0 
               tri_stress = pre_stress.' + d_o_vec;
               mean_stress = ( tri_stress(1)+tri_stress(2)+tri_stress(3) ) / 3;
               tri_dev_stress = [tri_stress(1) - mean_stress
                                 tri_stress(2) - mean_stress
                                 tri_stress(3) - mean_stress
                                 tri_stress(4)
                                 tri_stress(5)
                                 tri_stress(6)];
               tri_mises_stress = sqrt( 3/2 * (tri_dev_stress(1)^2 + tri_dev_stress(2)^2 + tri_dev_stress(3)^2 + ...
                                               2*tri_dev_stress(4)^2 + 2*tri_dev_stress(5)^2 + 2*tri_dev_stress(6)^2) );
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
                   final_stress = final_dev_stress + [mean_stress mean_stress mean_stress 0 0 0].';
                   %final_stress.'
                   %d_o_vec_tri = d_o_vec 
                   d_o_vec = final_stress - pre_stress.';
               end

          end

          d_integ_stress_(index_i,:) = d_o_vec.';
		  d_integ_strain_(index_i,:) = d_e_vec.';
       
          o_vec = pre_stress.' + d_o_vec;
          q_vec_i = B.' * o_vec;
          q_vec = q_vec + W(i)*W(i)*W(i) * detJarray(i) * q_vec_i;
        end %i = 1 : integ_num

        %Q_element_(e, :) = q_vec.';
        for i = 1 : 8
            %Q( [1:3] + (elementmat(e,i)-1)*3 ) = Q( [1:3] + (elementmat(e,i)-1)*3 ) + q_vec([1:3] + (i-1)*3 );
            Q( 1 + (elementmat(e,i)-1)*3 ) = Q( 1 + (elementmat(e,i)-1)*3 ) + q_vec(1 + (i-1)*3 );
            Q( 2 + (elementmat(e,i)-1)*3 ) = Q( 2 + (elementmat(e,i)-1)*3 ) + q_vec(2 + (i-1)*3 );
            Q( 3 + (elementmat(e,i)-1)*3 ) = Q( 3 + (elementmat(e,i)-1)*3 ) + q_vec(3 + (i-1)*3 );
        end

     end

  end


  function [BV_i, BVbar_i] = cal_BVbar(Pusai_i, J_i, detJ_i)

    BV_i = zeros(6,24);
    BVbar_i = zeros(6,24);

    %P2 = inv(J_i) * Pusai_i;
    P2 = J_i \ Pusai_i;
    N = reshape(P2,1,24);

    BV_i(1,:) = N;
    BV_i(2,:) = N;
    BV_i(3,:) = N;
    BV_i = BV_i / 3;
    %BV_i = [N; N; N; zeros(1,24); zeros(1,24); zeros(1,24)] / 3;

    BVbar_i = BV_i * detJ_i;

  end


  function [B, J] = cal_B_hexa(Pusai1, e_position )

    J = Pusai1 * e_position;
    %P2 = inv(J) * Pusai1;
    P2 = J \ Pusai1;

    B = zeros(6,24);
    % B(1,1:3:24) = P2(1,:);
    % B(2,2:3:24) = P2(2,:);
    % B(3,3:3:24) = P2(3,:);
    % B(4,1:3:24) = P2(2,:);
    % B(4,2:3:24) = P2(1,:);
    % B(5,2:3:24) = P2(3,:);
    % B(5,3:3:24) = P2(2,:);
    % B(6,1:3:24) = P2(3,:);
    % B(6,3:3:24) = P2(1,:); %-> slow

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

        nElement = size(elementmat,1);
        element_ctr = zeros(nElement,3);

        for i = 1 : nElement
            element_ctr(i,:) = sum(position(elementmat(i,:),:)) / 8;
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
            elem = MODEL.PART( part_id ).elementmat(j,:);
            faces(6*(c-1)+1, :) = [elem(1:4)];
            faces(6*(c-1)+2, :) = [elem(5:8)];
            faces(6*(c-1)+3, :) = [elem(1) elem(2) elem(6) elem(5)];
            faces(6*(c-1)+4, :) = [elem(2) elem(3) elem(7) elem(6)];
            faces(6*(c-1)+5, :) = [elem(3) elem(4) elem(8) elem(7)];
            faces(6*(c-1)+6, :) = [elem(4) elem(1) elem(5) elem(8)];
            faces_eleid(6*(c-1)+(1:6)) = j;
            ctr = sum( cdmat(elem,:) ) / 8;
            for k = 1 : 6
                index = 6*(c-1)+k;
                v1 = cdmat( faces(index, 2), :) - cdmat( faces(index, 1), :);
                v2 = cdmat( faces(index, 4), :) - cdmat( faces(index, 1), :);
                nv = cross(v1,v2);
                vc = ctr - cdmat( faces(index, 1), :);
                if nv * vc.' > 0
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


  function [c_triangles, c_triangles_eleid, c_nodes] = get_surface_triangle(MODEL, i, array_element)

        nE = length(array_element); %MODEL.PART( part_id ).nElement;
        surfaces = zeros(nE*6, 4);
        surfaces_eleid = zeros(nE*6, 1);
        sorted_surfaces = zeros(nE*6, 4);
        c = 1;
        for j = array_element
            for k = 1 : 6
                surfaces(6*(c-1)+k, :) = MODEL.INSTANCE(i).surfaces(6*(j-1)+k,:);
                sorted_surfaces(6*(c-1)+k, :) = MODEL.INSTANCE(i).sorted_surfaces(6*(j-1)+k,:);
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
        %MODEL.INSTANCE(i).c_triangles = c_triangles;
        %MODEL.INSTANCE(i).c_triangles_eleid = c_triangles_eleid;

        c_nodes = reshape(c_surfaces, 1, size(c_surfaces,1)*4);
        c_nodes = unique(c_nodes);
        %length(c_nodes);
        %MODEL.INSTANCE(i).c_nodes = c_nodes;

  end

  function c_force3 = cal_contact_force(c_pair, position, velo, diag_M, elementMinSize, ...
                                        INSTANCE, PART, MATERIAL, element_flag)

    nNode = size(position,1);
    c_force3 = zeros(nNode,3);   

    %'cal_contact_force'

    for c = 1 : size(c_pair,1)
        i_part = c_pair(c,1);
        j_part = c_pair(c,2);
        i_part_id = INSTANCE(i_part).part_id;
        j_part_id = INSTANCE(j_part).part_id;

        node_offset_i = INSTANCE(i_part).node_offset;
        node_offset_j = INSTANCE(j_part).node_offset;
        nNode_i = INSTANCE(i_part).nNode;
        nNode_j = INSTANCE(j_part).nNode;
        c_nodes = INSTANCE(i_part).c_nodes + node_offset_i;
        element_offset_j = INSTANCE(j_part).element_offset;

        c_triangles = INSTANCE(j_part).c_triangles + node_offset_j;
        c_triangles_eleid = INSTANCE(j_part).c_triangles_eleid + element_offset_j;
        young = MATERIAL( PART(j_part_id).material_id ).young;
        d_min = elementMinSize * 0.3;

        %--- contact range ---%
        position_i = position((1:nNode_i)+node_offset_i,:);
        position_j = position((1:nNode_j)+node_offset_j,:);
        min_i = min(position_i);
        max_i = max(position_i);
        min_j = min(position_j);
        max_j = max(position_j);

        range_min_x = max(min_i(1),min_j(1));
        range_max_x = min(max_i(1),max_j(1));
        range_min_y = max(min_i(2),min_j(2));
        range_max_y = min(max_i(2),max_j(2));
        range_min_z = max(min_i(3),min_j(3));
        range_max_z = min(max_i(3),max_j(3));

        if range_min_x > range_max_x || range_min_y > range_max_y || range_min_z > range_max_z
            continue
        end


        for j = 1 : size(c_triangles,1)
            if element_flag( c_triangles_eleid(j) ) == 0
                continue
            end
            j0 = c_triangles(j,1);
            j1 = c_triangles(j,2);
            j2 = c_triangles(j,3);
            q0 = position(j0,:);
            q1 = position(j1,:);
            q2 = position(j2,:);

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
            

            v1 = q1 -q0;
            v2 = q2 -q0;
            L1 = norm(v1);
            L2 = norm(v2);
            Lmax = max(L1,L2);
            n = cross(v1,v2);
            n = n / norm(n);
            S = 0.5 * sqrt( v1*v1.' * v2*v2.' - (v1*v2.')^2 );
            L = (norm(v1) + norm(v2))*0.5 * 2.0;

            %--- custom
            % if abs(n(3)) > 0.99
            %     continue
            % end

            A = zeros(3,3);
            A(:,1) = v1.';
            A(:,2) = v2.';
            A(:,3) = -n.';

            for i = c_nodes
                p = position(i,:);

                if p(1) < range_min_x || p(2) < range_min_y || p(3) < range_min_z
                    continue
                end
                if p(1) > range_max_x || p(2) > range_max_y || p(3) > range_max_z
                    continue
                end


                b = (p - q0).';
                if norm(b) > L
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
                
                x = A \ b;
                
                if 0.0 <= x(1) && 0.0 <= x(2) && x(1) + x(2) <= 1.0 && ...
                        x(3) > 0 && x(3) <= d_min
                     % x
                     % i
                     % p
                     % q0
                     % q1
                     % q2
                     % n
                     % dA
                     %x(3)
                    
                    %v = [velo(i*3-2) velo(i*3-1) velo(i*3)];
                    k = young * S / Lmax;
                    %cc = 2 * sqrt( diag_M(i) * k ) * 0;
                    f = k * x(3) * n; %  - cc * (v*n.') * n;
                    c_force3(i,:) = c_force3(i,:) + f; 
                    c_force3(j0,:) = c_force3(j0,:) - f  / 3;
                    c_force3(j1,:) = c_force3(j1,:) - f  / 3;
                    c_force3(j2,:) = c_force3(j2,:) - f  / 3;
                    % if sum(abs(imag(f))) > 0
                    %  x
                    %  i
                    %  p
                    %  q0
                    %  q1
                    %  q2
                    %  v1
                    %  v2
                    %  S
                    %  n
                    %  det(A)
                    % end
                end

            end
        end

    end

  end


  function [node_value] = cal_node_stress_strain(nNode, elementmat, integ_num, data)

    integ_stress = data.integ_stress;
    integ_strain = data.integ_strain;
    integ_eq_plastic_strain = data.integ_eq_plastic_strain;
    integ_triax_stress = data.integ_triax_stress;

    node_stress = zeros(nNode, 6);
    node_strain = zeros(nNode, 6);
    node_eq_plastic_strain = zeros(nNode,1);
    node_triax_stress = zeros(nNode,1);

    nElement = size(elementmat,1);
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
            element_stress(i,:) = sum( integ_stress([1:integ_num]+(i-1)*integ_num, :)  ) / integ_num;
            element_strain(i,:) = sum( integ_strain([1:integ_num]+(i-1)*integ_num, :)  ) / integ_num;
            element_eq_plastic_strain(i) = sum( integ_eq_plastic_strain([1:integ_num]+(i-1)*integ_num)  ) / integ_num;
            element_triax_stress(i) = sum( integ_triax_stress([1:integ_num]+(i-1)*integ_num)  ) / integ_num;
        end
    end

    for i = 1 : nElement
        for k = 1 : 8
            node_stress(elementmat(i,k),:) = node_stress(elementmat(i,k),:) + element_stress(i, :);
            node_strain(elementmat(i,k),:) = node_strain(elementmat(i,k),:) + element_strain(i, :);
            node_eq_plastic_strain(elementmat(i,k)) = node_eq_plastic_strain(elementmat(i,k)) + element_eq_plastic_strain(i);
            node_triax_stress(elementmat(i,k)) = node_triax_stress(elementmat(i,k)) + element_triax_stress(i);
        end
    end

    inc_num = zeros(nNode, 1);
    for i = 1 : nElement
       inc_num(elementmat(i,:)) = inc_num(elementmat(i,:)) + 1;
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


function write_vtk(index_, coordmat, elementmat, element_flag, disp, data)

    node_stress = data.node_stress;
    node_strain = data.node_strain;
    node_mises_stress = data.node_mises_stress;
    node_eq_plastic_strain = data.node_eq_plastic_strain;
    node_triax_stress = data.node_triax_stress;

    nNode = size(coordmat,1);
    nElement = size(elementmat,1);
    disp3 = reshape(disp,3,nNode)';

    for i = 1 : nNode
        for j = 1 : 3
            if abs(disp3(i,j)) < 1E-16
                disp3(i,j) = 0;
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

    fname = sprintf("temp\\file%03d.vtk", index_);
    out = fopen(fname,"w");
    fprintf(out,"# vtk DataFile Version 2.0\n");
    fprintf(out,"Test\n");
    fprintf(out,"ASCII\n");
    fprintf(out,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(out,"POINTS %d float\n", nNode);
    
    for i = 1 : nNode
        fprintf(out, sprintf("%1.6e %1.6e %1.6e\n", coordmat(i,1), coordmat(i,2), coordmat(i,3)) );
    end

    draw_element_num = sum(element_flag);

    %fprintf(out, "CELLS %d %d\n", nElement, nElement*(8+1));
    fprintf(out, "CELLS %d %d\n", draw_element_num, draw_element_num*(8+1));
    e = elementmat;
    e = e - 1;
    for i = 1 : nElement
        if element_flag(i) == 1
            fprintf(out, sprintf("8 %d %d %d %d %d %d %d %d\n", ...
            e(i,1), e(i,2), e(i,3), e(i,4), e(i,5), e(i,6), e(i,7), e(i,8)) );
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

