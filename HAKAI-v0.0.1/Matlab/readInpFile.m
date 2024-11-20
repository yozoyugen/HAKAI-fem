function [MODEL] = readInpFile(fname)

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


  fid = fopen(fname);

  lines = {};
  n = 1;
    while ~feof(fid)
        tline = fgetl(fid);
        lines{n} = tline;
        n = n + 1;
    end
  n = n-1;

  fclose(fid);
  %lines

  
  %--- Part ---%
  part_index = [];
  nPart = 0;
  for i = 1 : n
    if strfind(lines{i}, '*Part, name=') > 0
        part_index = [part_index i];
        nPart = nPart + 1;
    end
  end
  PART(1:nPart) = struct('name', "");

  for k = 1 : nPart
      s = replace(lines{part_index(k)},' ', '');
      ss = split( s,',');
      ss = ss{2};
      name_ = ss( (strfind(ss, 'name=') + 5):length(ss) );
      PART(k).name = name_;

      %--- Node ---%
      index = 1;
      for i = part_index(k) : n
        if strfind(lines{i}, '*Node') > 0
            index = i;
            break;
        end
      end

      nNode = 0;
      for i = index+1 : n
        if strfind(lines{i}, '*') > 0
            break;
        end
        nNode = nNode + 1;
      end
      nNode
      PART(k).nNode = nNode;
      coordmat = zeros(nNode, 3);

      for i = 1 : nNode
        if strfind(lines{index+i}, '*') > 0
            break;
        end
        s = replace(lines{index+i},' ', '');
        ss = split(s,',');
        coordmat(i,1) = str2double( ss{2});
        coordmat(i,2) = str2double( ss{3});
        coordmat(i,3) = str2double( ss{4});
      end
      PART(k).coordmat = coordmat'; % column major

      %--- Element ---%
      index = 1;
      for i = part_index(k) : n
        if strfind(lines{i}, '*Element') > 0
            index = i;
            s = replace(lines{i},' ', '');
            ss = split( s,',');
            ss = ss{2};
            name_ = ss( (strfind(ss, 'type=') + 5):length(ss) );
            PART(k).element_type = name_;
            break;
        end
      end
    
      nElement = 0;
      for i = index+1 : n
        if strfind(lines{i}, '*') > 0
            break;
        end
        nElement = nElement + 1;
      end
      nElement
      PART(k).nElement = nElement;
    
      jn = 8;
      if contains(PART(k).element_type, 'S4') > 0  
          jn = 4;
      end
      elementmat = zeros(nElement, jn);

      for i = 1 : nElement
        if strfind(lines{index+i}, '*') > 0
            break;
        end
        s = replace(lines{index+i},' ', '');
        ss = split(s,',');
        for j = 1 : jn
            elementmat(i,j) = str2num( ss{1+j});
        end
      end
      %elementmat
      PART(k).elementmat = elementmat'; %column major

      %--- Nset ---%
      nset_index = [];
      for i = part_index(k) : n
        if contains(lines{i}, '*Nset') > 0  && contains(lines{i}, 'generate') > 0
            nset_index = [nset_index i];
        end
        if contains(lines{i}, '*End Part') > 0
            break
        end
      end
    
      nset_num = length(nset_index)
      NSET_p(1:nset_num) = struct('name', "");
    
      for i = 1 : nset_num
          s = replace(lines{nset_index(i)},' ', '');
          ss = split( s,',');
          ss = ss{2};
          nset_name = ss( (strfind(ss, 'nset=') + 5):length(ss) );
          NSET_p(i).name = nset_name;
          
          s = replace(lines{nset_index(i)+1},' ', '');
          ss = split( s,',');
          a = [str2num(ss{1}) : str2num(ss{3}) : str2num(ss{2})];
          NSET_p(i).nodes = a;
      end
      PART(k).NSET = NSET_p;
      clear NSET_p;

      for i = part_index(k) : n
        if contains(lines{i}, '*Solid Section') > 0  || contains(lines{i}, '*Shell Section')
          s = replace(lines{i},' ', '');
          ss = split( s,',');
          for j = 1 : length(ss)
              sss = ss{j};
              if strfind(sss, 'material=') > 0
                name_ = sss( (strfind(sss, 'material=') + 9):length(sss) );
                PART(k).material_name = name_;
                break
              end
          end

          if contains(lines{i}, '*Shell') 
              s = replace(lines{i+1},' ', '');
              ss = split( s,',');
              PART(k).shell_thickness = str2num(ss{1});
          end

          break
        end
      end
  end


  %--- Instance ---%
  instance_index = [];
  for i = 1 : n
    if contains(lines{i}, '*Instance') > 0  
        instance_index = [instance_index i];
    end
  end
  instance_num = length(instance_index);
  INSTANCE(1:instance_num) = struct('name', "");

  for k = 1 : instance_num
      s = replace(lines{instance_index(k)},' ', '');
      ss = split( s,',');
      ss = ss{2};
      name_ = ss( (strfind(ss, 'name=') + 5):length(ss) );
      INSTANCE(k).name = name_;

      ss = split( s,',');
      ss = ss{3};
      name_ = ss( (strfind(ss, 'part=') + 5):length(ss) );
      INSTANCE(k).part_name = name_;

      for i = 1 : nPart
          if strcmp(PART(i).name, INSTANCE(k).part_name )
            INSTANCE(k).part_id = i;
            break
          end
      end

      INSTANCE(k).translate = {};
      j = 1;
      for i = instance_index(k)+1 : n
        if contains(lines{i}, '*End Instance') > 0  
            break
        end

          s = replace(lines{i},' ', '');
          INSTANCE(k).translate{j} = s;
          j = j + 1;
      end

  end
  
  %--- Nset ---%
  nset_index = [];
  for i = 1 : n
    if contains(lines{i}, '*Nset') > 0  && contains(lines{i}, 'instance=') > 0
        nset_index = [nset_index i];
    end
  end

  nset_num = length(nset_index);
  %NSET = {};
  NSET(1:nset_num) = struct('name', "");

  for k = 1 : nset_num
      s = replace(lines{nset_index(k)},' ', '');
      ss = split( s,',');
      ss = ss{2};
      nset_name = ss( (strfind(ss, 'nset=') + 5):length(ss) );
      NSET(k).name = nset_name;
     
      ss = split( s,',');
      ss = ss{3};
      name_ = ss( (strfind(ss, 'instance=') + 9):length(ss) );
      NSET(k).instance_name = name_;

      for i = 1 : instance_num
          if strcmp(NSET(k).instance_name, INSTANCE(i).name)
                NSET(k).part_name = INSTANCE(i).part_name;
                NSET(k).part_id = INSTANCE(i).part_id;
                NSET(k).instance_id = i;
          end
      end

      a = [];
      index = nset_index(k);

      ss = split( s,',');
      if length(ss) == 4 && strcmp(ss{4},"generate")

           s = replace(lines{index+1},' ', '');
           ss = split( s,',');
           a = [str2num(ss{1}) : str2num(ss{3}) : str2num(ss{2})];
           
      else

          for i = index+1 : n
            if strfind(lines{i}, '*') > 0
                break;
            end
    
            s = replace(lines{i},' ', '');
            ss = split( s,',');
            for j = 1 : length(ss)
                if length(ss{j}) > 0
                    a = [a  str2num( ss{j})];
                end
            end
          end

      end
      NSET(k).nodes = a;
  end


  %--- Elset ---%
  elset_index = [];
  for i = 1 : n
    if contains(lines{i}, '*Elset') > 0  && contains(lines{i}, 'instance=') > 0
        elset_index = [elset_index i];
    end
  end

  elset_num = length(elset_index);
  ELSET(1:elset_num) = struct('name', "");

  for k = 1 : elset_num
      s = replace(lines{elset_index(k)},' ', '');
      ss = split( s,',');
      ss = ss{2};
      name_ = ss( (strfind(ss, 'elset=') + 6):length(ss) );
      ELSET(k).name = name_;
     
      ss = split( s,',');
      ss = ss{3};
      if contains(ss, 'instance=') > 0
          name_ = ss( (strfind(ss, 'instance=') + 9):length(ss) );
          ELSET(k).instance_name = name_;
      else
          ss = split( s,',');
          ss = ss{4};
          if contains(ss, 'instance=') > 0
              name_ = ss( (strfind(ss, 'instance=') + 9):length(ss) );
              ELSET(k).instance_name = name_;
          end
      end

      for i = 1 : instance_num
          if strcmp(ELSET(k).instance_name, INSTANCE(i).name)
                ELSET(k).part_name = INSTANCE(i).part_name;
                ELSET(k).part_id = INSTANCE(i).part_id;
                ELSET(k).instance_id = i;
          end
      end

      a = [];
      index = elset_index(k);

      ss = split( s,',');
      if length(ss) == 4 && strcmp(ss{4},"generate")

           s = replace(lines{index+1},' ', '');
           ss = split( s,',');
           a = str2num(ss{1}) : str2num(ss{3}) : str2num(ss{2});

      elseif length(ss) == 5 && strcmp(ss{3},"internal") && strcmp(ss{5},"generate")

           s = replace(lines{index+1},' ', '');
           ss = split( s,',');
           a = str2num(ss{1}) : str2num(ss{3}) : str2num(ss{2});

      elseif length(ss) == 4 && strcmp(ss{3},"internal")

          for i = index+1 : n
            if strfind(lines{i}, '*') > 0
                break;
            end

            s = replace(lines{i},' ', '');
            ss = split( s,',');
            for j = 1 : length(ss)
                if length(ss{j}) > 0
                    a = [a  str2num( ss{j})];
                end
            end
          end

      end
      ELSET(k).elements = a;

  end


  %--- Surface ---%
  surface_index = [];
  for i = 1 : n
    if contains(lines{i}, '*Surface,') > 0  
        surface_index = [surface_index i];
    end
  end

  surface_num = length(surface_index);
  SURFACE(1:surface_num) = struct('name', "");

  for k = 1 : surface_num
      s = replace(lines{surface_index(k)},' ', '');
      ss = split( s,',');
      ss = ss{3};
      name_ = ss( (strfind(ss, 'name=') + 5):length(ss) );
      SURFACE(k).name = name_;

      index = surface_index(k);
      a = [];
      c = 1;
      for i = index+1 : n
        if strfind(lines{i}, '*') > 0
            break;
        end

        s = replace(lines{i},' ', '');
        ss = split( s,',');
        SURFACE(k).elset_name{c} = ss{1};
        
        for j = 1 : elset_num
            if strcmp(SURFACE(k).elset_name{c}, ELSET(j).name)
                SURFACE(k).instance_id = ELSET(j).instance_id;
                a = [a ELSET(j).elements];
            end
        end

        c = c + 1;
      end
      SURFACE(k).elements = unique(sort(a));
         
  end
  
  

  %--- Global Model ---%
  nNode = 0;
  nElement = 0;
  coordmat = [];
  elementmat = [];
  for i = 1 : instance_num
      part_id = INSTANCE(i).part_id;
      coordmat_i = PART(part_id).coordmat;
      INSTANCE(i).node_offset = nNode;
      INSTANCE(i).element_offset = nElement;
      INSTANCE(i).nNode = PART(part_id).nNode;
      INSTANCE(i).nElement = PART(part_id).nElement;
      

      T =  eye(3);
      for j = length(INSTANCE(i).translate): -1 : 1
          s = INSTANCE(i).translate{j};
          ss = split( s,',');
          if length(ss) == 3
              offset_ = [str2num(ss{1}) str2num(ss{2}) str2num(ss{3})]';
              coordmat_i = coordmat_i + offset_ * ones(1,size(coordmat_i,2));
          elseif length(ss) == 7
              nv = [str2num(ss{4})-str2num(ss{1}) 
                    str2num(ss{5})-str2num(ss{2}) 
                    str2num(ss{6})-str2num(ss{3}) ];
              nv = nv / norm(nv);
              n1 = nv(1);
              n2 = nv(2);
              n3 = nv(3);
              d = str2num(ss{7}) / 180.0 * pi;
              T = [n1*n1*(1-cos(d))+cos(d)  n1*n2*(1-cos(d))-n3*sin(d)  n1*n3*(1-cos(d))+n2*sin(d)
                   n1*n2*(1-cos(d))+n3*sin(d)  n2*n2*(1-cos(d))+cos(d)  n2*n3*(1-cos(d))-n1*sin(d)
                   n1*n3*(1-cos(d))-n2*sin(d)  n2*n3*(1-cos(d))+n1*sin(d)  n3*n3*(1-cos(d))+cos(d) ];
              coordmat_i = T * coordmat_i;
          end
      end

      coordmat = [coordmat  coordmat_i];  % column major
      
      elementmat_i = PART(part_id).elementmat + nNode;
      
      if size(elementmat_i,1) == 4
            elementmat_i = [elementmat_i; zeros(size(elementmat_i))];
      end

      elementmat = [elementmat  elementmat_i]; % column major
      
      nNode = nNode + PART(part_id).nNode;
      nElement = nElement + PART(part_id).nElement;
  end


  %--- Amplitude ---%
  amplitude_index = [];
  for i = 1 : n
    if contains(lines{i}, '*Amplitude') > 0  
        amplitude_index = [amplitude_index i];
    end
  end

  amplitude_num = length(amplitude_index);
  AMPLITUDE(1:amplitude_num) = struct('name', "", 'time', [], 'value', []);

  for k = 1 : amplitude_num
      s = replace(lines{amplitude_index(k)},' ', '');
      ss = split( s,',');
      ss = ss{2};
      name_ = ss( (strfind(ss, 'name=') + 5):length(ss) );
      AMPLITUDE(k).name = name_;
      
      index = amplitude_index(k);  
      for i = index+1 : n
        if strfind(lines{i}, '*') > 0
            break;
        end
        
        s = replace(lines{i},' ', '');
        ss = split( s,',');
        for j = 1 : length(ss)/2
            AMPLITUDE(k).time(j) = str2num(ss{j*2-1});
            AMPLITUDE(k).value(j) = str2num(ss{j*2});
        end

      end

  end


  %--- Material ---%
  material_index = [];
  for i = 1 : n
    if contains(lines{i}, '*Material') > 0  
        material_index = [material_index i];
    end
  end

  material_num = length(material_index);
  METERIAL(1:material_num) = struct('name', "");

  for k = 1 : material_num
      s = replace(lines{material_index(k)},' ', '');
      ss = split( s,',');
      ss = ss{2};
      name_ = ss( (strfind(ss, 'name=') + 5):length(ss) );
      MATERIAL(k).name = name_;
      MATERIAL(k).fracture_flag = 0;
      MATERIAL(k).failure_stress = [];

      index = material_index(k);
      plastic_index = 1;
      ductile_index = 1;

      for i = index+1 : n
        if strfind(lines{i}, '*Material') > 0
            break;
        end
        if strfind(lines{i}, '**') > 0
            break;
        end

        if strfind(lines{i}, '*Density') > 0
            s = replace(lines{i+1},' ', '');
            ss = split( s,',');
            MATERIAL(k).density = str2num( ss{1} );
        end

        if strfind(lines{i}, '*Elastic') > 0
            s = replace(lines{i+1},' ', '');
            ss = split( s,',');
            MATERIAL(k).young = str2num( ss{1} );
            MATERIAL(k).poisson = str2num( ss{2} );
        end

        if strfind(lines{i}, '*Plastic') > 0
            plastic_index = i;
        end

        if contains(lines{i}, '*Damage Initiation')  && contains(lines{i}, 'criterion=DUCTILE')
            ductile_index = i
            MATERIAL(k).fracture_flag = 1;
        end
        
        if strfind(lines{i}, '*Tensile Failure') > 0
            s = replace(lines{i+1},' ', '');
            ss = split( s,',');
            MATERIAL(k).failure_stress = str2num( ss{1} );
            MATERIAL(k).fracture_flag = 1;
        end

      end

      plastic = [];
      if plastic_index > material_index(k)
          for i = plastic_index+1 : n
            if strfind(lines{i}, '*') > 0
                break;
            end
            s = replace(lines{i},' ', '');
            ss = split( s,',');
            plastic = [plastic; str2num(ss{1})  str2num(ss{2})];
          end
      end
      MATERIAL(k).plastic = plastic;


      ductile = [];
      if ductile_index > material_index(k)
          for i = ductile_index+1 : n
            if strfind(lines{i}, '*') > 0
                break;
            end
            s = replace(lines{i},' ', '');
            ss = split( s,',');
            ductile = [ductile; str2num(ss{1})  str2num(ss{2})  str2num(ss{3})];
          end
      end
      MATERIAL(k).ductile = ductile;
      
  end

  element_material = [];
  element_instance = [];
  for i = 1 : length(INSTANCE) %length(PART)
      part_id = INSTANCE(i).part_id;
       for j = 1 : length(MATERIAL)
           % if strcmp(PART(i).material_name, MATERIAL(j).name)
           %      PART(i).material_id = j;
           % end
           if strcmp(PART(part_id).material_name, MATERIAL(j).name)
                PART(part_id).material_id = j;
                INSTANCE(i).material_id = j;
           end
       end
       mat_id = ones(PART(part_id).nElement,1) * PART(part_id).material_id;
       element_material = [element_material; mat_id];

       instance_id = ones(PART(part_id).nElement,1) * i;
       element_instance = [element_instance; instance_id];
  end


  %--- Step ---%
  d_time = 0;
  end_time = 0;
  for i = 1 : n
    if contains(lines{i}, '*Dynamic, Explicit') > 0  
        s = replace(lines{i+1},' ', '');
        ss = split( s,',');
        d_time = str2double(ss{1});
        end_time = str2double(ss{2});
        break
    end
  end

  %--- Mass scaling ---%
  mass_scaling = 1;
  for i = 1 : n
    if contains(lines{i}, '*Fixed Mass Scaling') > 0  
        s = replace(lines{i},' ', '');
        ss = split( s,',');
        ss = ss{2};
        v_ = ss( (strfind(ss, 'factor=') + 7):length(ss) );
        mass_scaling = str2double(v_);
        break
    end
  end


  %--- BC ---%
  bc_index = [];
  for i = 1 : n
    if contains(lines{i}, '*Boundary') > 0  
        bc_index = [bc_index i];
    end
  end

  bc_num = length(bc_index);
  BC(1:bc_num) = struct('Nset_name', '',  'dof', [], 'value', [], ...
                        'amp_name', '', 'amplitude', []);

  for k = 1 : bc_num
      index = bc_index(k);
      %BC(k).dof = [];

      s = replace(lines{index},' ', '');
      ss = split( s,',');
      if length(ss) == 2 && contains(ss{2}, 'amplitude=')
          ss = ss{2};
          name_ = ss( (strfind(ss, 'amplitude=') + 10):length(ss) );
          BC(k).amp_name = name_;
          for j = 1 : amplitude_num
                if strcmp(AMPLITUDE(j).name, name_)
                    BC(k).amplitude = AMPLITUDE(j);
                    break
                end
          end
      end


      %dof = [];
      for i = index+1 : n
        if strfind(lines{i}, '*Boundary') > 0
            break;
        end
        if strfind(lines{i}, '**') > 0
            break;
        end

        s = replace(lines{i},' ', '');
        ss = split( s,',');
        BC(k).Nset_name = ss{1};

        nodes = [];
        if strfind(BC(k).Nset_name, '.') > 0
            sss = split( BC(k).Nset_name,'.');
            instance_name = sss{1};
            nset_name = sss{2};
            instance_id = 0;
            part_id = 0;
            for j = 1 : instance_num
                if strcmp(INSTANCE(j).name, instance_name)
                    instance_id = j;
                    part_id = INSTANCE(j).part_id;
                    break
                end
            end

            for j = 1 : length(PART(part_id).NSET)
                if strcmp(PART(part_id).NSET(j).name, nset_name)
                    nodes_j = PART(part_id).NSET(j).nodes + INSTANCE(instance_id).node_offset;
                    nodes = [nodes  nodes_j];
                    break
                end
            end

        else 
            for j = 1 : nset_num
                if strcmp(BC(k).Nset_name, NSET(j).name)
                    nodes_j = NSET(j).nodes + INSTANCE(NSET(j).instance_id).node_offset;
                    nodes = [nodes nodes_j];
                end
            end
        end

        dof = [];
        if length(ss) == 2 && strcmp(ss{2}, 'ENCASTRE')
            dof = [dof nodes*3-2 nodes*3-1 nodes*3-0 ];
            BC(k).dof{1} = sort(dof);
            BC(k).value{1} = 0;
        elseif length(ss) == 3
            ii = str2num(ss{2});
            dir = str2num(ss{3});
            if dir <= 3
                dof = [nodes*3-(3-dir) ];
                BC(k).dof{ii} = dof;
                BC(k).value{ii} = 0;
            end
        elseif length(ss) == 4
            ii = str2num(ss{2});
            dir = str2num(ss{3});
            value = str2num(ss{4});
            if dir <= 3
                dof = [dof nodes*3-(3-dir) ];
                BC(k).dof{ii} = dof;
                BC(k).value{ii} = value;
            end

        end

        %BC(k).dof = [BC(k).dof  sort(dof)];

      end

      

  end


  %--- Initial Conditions ---%
  ic_index = [];
  for i = 1 : n
    if contains(lines{i}, '*Initial Conditions') > 0  
        ic_index = [ic_index i];
    end
  end

  ic_num = length(ic_index);
  IC(1:ic_num) = struct('Nset_name', '',  'dof', []);

  for k = 1 : ic_num
      index = ic_index(k);

      s = replace(lines{index},' ', '');
      ss = split( s,',');
      ss = ss{2};
      type_ = ss( (strfind(ss, 'type=') + 5):length(ss) );
      IC(k).type = type_;
      
      for i = index+1 : n
        if strfind(lines{i}, '*Initial Conditions') > 0
            break;
        end
        if strfind(lines{i}, '**') > 0
            break;
        end

        s = replace(lines{i},' ', '');
        ss = split( s,',');
        IC(k).Nset_name = ss{1};

        nodes = [];
        if strfind(IC(k).Nset_name, '.') > 0
            sss = split( IC(k).Nset_name,'.');
            instance_name = sss{1};
            nset_name = sss{2};
            instance_id = 0;
            part_id = 0;
            for j = 1 : instance_num
                if strcmp(INSTANCE(j).name, instance_name)
                    instance_id = j;
                    %part_id = j;
                    part_id = INSTANCE(j).part_id;
                    break
                end
            end

            for j = 1 : length(PART(part_id).NSET)
                if strcmp(PART(part_id).NSET(j).name, nset_name)
                    nodes_j = PART(part_id).NSET(j).nodes + INSTANCE(instance_id).node_offset;
                    nodes = [nodes  nodes_j];
                    break
                end
            end

        else 
            for j = 1 : nset_num
                if strcmp(IC(k).Nset_name, NSET(j).name)
                    nodes_j = NSET(j).nodes + INSTANCE(NSET(j).instance_id).node_offset;
                    nodes = [nodes nodes_j];
                    break
                end
            end
        end

        dof = [];
        dir = str2num(ss{2});
        dof = [dof nodes*3-(3-dir) ];

        v = str2num(ss{3});
        
        IC(k).dof = [IC(k).dof  dof];
        IC(k).value = v;

      end

  end

  %--- Interactions (Contact) ---%
  contact_flag = 0;
  for i = 1 : n
    if contains(lines{i}, '*Contact') > 0   %'*Contact Inclusions'
        contact_flag = 1;
        break
    end
  end

  for i = 1 : n
    if contains(lines{i}, '*Contact Inclusions') > 0  && contains(lines{i}, 'HAKAIoption=self-contact') > 0
        contact_flag = 2;
        break
    end
  end

  cp_index = [];
  for i = 1 : n
    if contains(lines{i}, '*Contact Pair,') > 0  
        cp_index = [cp_index i];
    end
  end

  cp_num = length(cp_index);
  CP(1:cp_num) = struct('name', '');
  
  for k = 1 : cp_num
      index = cp_index(k);
      s = replace(lines{index},' ', '');
      ss = split( s,',');
      ss = ss{4};
      name_ = ss( (strfind(ss, 'cpset=') + 6):length(ss) );
      CP(k).name = name_;

      s = replace(lines{index+1},' ', '');
      ss = split( s,',');
      CP(k).surface_name_i = ss{1};
      CP(k).surface_name_j = ss{2};

        for j = 1 : surface_num
            if strcmp(CP(k).surface_name_i, SURFACE(j).name)
                CP(k).instance_id_i =  SURFACE(j).instance_id;
                CP(k).elements_i = SURFACE(j).elements;
            end
            if strcmp(CP(k).surface_name_j, SURFACE(j).name)
                CP(k).instance_id_j =  SURFACE(j).instance_id;
                CP(k).elements_j = SURFACE(j).elements;
            end
        end

  end


  MODEL.PART = PART;
  MODEL.INSTANCE = INSTANCE;
  MODEL.NSET = NSET;
  MODEL.ELSET = ELSET;
  MODEL.SURFACE = SURFACE;
  MODEL.AMPLITUDE = AMPLITUDE;
  MODEL.MATERIAL = MATERIAL;
  MODEL.BC = BC;
  MODEL.IC = IC;
  MODEL.CP = CP;
  MODEL.nNode = nNode;
  MODEL.coordmat = coordmat;
  MODEL.nElement = nElement;
  MODEL.elementmat = elementmat;
  MODEL.element_material = element_material;
  MODEL.element_instance = element_instance;
  MODEL.d_time = d_time;
  MODEL.end_time = end_time;
  MODEL.mass_scaling = mass_scaling;
  MODEL.contact_flag = contact_flag;

end
