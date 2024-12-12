#=
HAKAI - A 3-dimensional finite element program
Copyright (c) 2024 Yozo Yugen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.
=#



using LinearAlgebra

mutable struct NsetType
    name::String
    instance_name::String
    instance_id::Int
    part_name::String
    part_id::Int
    nodes::Array{Int,1}
end

mutable struct ELsetType
    name::String
    instance_name::String
    instance_id::Int
    part_name::String
    part_id::Int
    elements::Array{Int,1}
end

mutable struct SurfaceType
    name::String
    elset_name::Array{String,1}
    instance_id::Int
    elements::Array{Int,1}
end

mutable struct PartType
    name::String
    nNode::Int
    coordmat::Array{Float64,2}
    nElement::Int
    elementmat::Array{Int,2}
    NSET::Array{NsetType,1}
    material_name::String
    material_id::Int
end

mutable struct InstanceType
    name::String
    part_name::String
    part_id::Int
    material_id::Int
    translate::Array{String,1}
    node_offset::Int
    nNode::Int
    element_offset::Int
    nElement::Int
    elements::Array{Int,1}
    surfaces::Array{Int,2}
    sorted_surfaces::Array{Int,2}
    surfaces_eleid::Array{Int,1}
    c_triangles::Array{Int,2}
    c_triangles_eleid::Array{Int,1}
    c_nodes::Array{Int,1}
end

mutable struct AmplitudeType
    name::String
    time::Array{Float64,1}
    value::Array{Float64,1}
end

mutable struct MaterialType
    name::String
    density::Float64
    young::Float64
    poisson::Float64
    plastic::Array{Float64,2}  #Vector{Any} -> become very slow!!
    Hd::Array{Float64,1}
    fracture_flag::Int
    failure_stress::Float64
    ductile::Array{Float64,2}
    G::Float64
    Dmat::Array{Float64,2}
end

mutable struct BCType
    Nset_name::String
    dof::Array{Array{Int,1},1}    #Vector{Any}   
    value::Vector{Float64} 
    amp_name::String
    amplitude::AmplitudeType
end

mutable struct ICType
    Nset_name::String
    type::String
    dof::Array{Array{Int,1},1}   #Vector{Any} 
    value::Vector{Float64}
end

mutable struct CPType
    name::String
    surface_name_1::String
    surface_name_2::String
    instance_id_1::Int
    instance_id_2::Int
    elements_1::Array{Int,1}
    elements_2::Array{Int,1}
    c_triangles_1::Array{Int,2}
    c_triangles_2::Array{Int,2}
    c_triangles_eleid_1::Array{Int,1}
    c_triangles_eleid_2::Array{Int,1}
    c_nodes_1::Array{Int,1}
    c_nodes_2::Array{Int,1}
end

mutable struct ModelType
    PART::Array{PartType,1}
    INSTANCE::Array{InstanceType,1}
    NSET::Array{NsetType,1}
    ELSET::Array{ELsetType,1}
    SURFACE::Array{SurfaceType,1}
    AMPLITUDE::Array{AmplitudeType,1}
    MATERIAL::Array{MaterialType,1}
    BC::Array{BCType,1}
    IC::Array{ICType,1}
    CP::Array{CPType,1}
    nNode::Int
    coordmat::Array{Float64,2}
    nElement::Int
    elementmat::Array{Int,2}
    element_material::Vector{Int}  #Vector{Any} -> become very slow!!
    element_instance::Vector{Int}
    d_time::Float64
    end_time::Float64
    mass_scaling::Float64
    contact_flag::Int
end

function readInpFile( fname )

    println("readInpFile:", fname) 

    f = open( fname, "r")
    lines = readlines(f)
　  close(f)

    #println(list) 
    n = length(lines)
    #println("n:", n) 
    #println("lines[1]:", lines[1]) 
    
    #--- Part ---#
    part_index = [];
    nPart = 0;
    for i = 1 : n
        #println("lines[i]:", lines[i]) 
        if occursin.("*Part, name=", lines[i]) 
            #println("i:", i) 
            push!(part_index, i)
            nPart = nPart + 1;
        end
    end
    #println("nPart:", nPart) 
    
    PART = PartType[]
    #p = PartType("bone", 1)
    #p.name = "beam"
    #push!(PART,p)
    
    for k = 1 : nPart

        p = PartType("", 0, zeros(1,1), 0, zeros(Int,1,1), NsetType[], "", 0)
        push!(PART,p)

        #println("part_index[k]:", part_index[k]) 
        #println("lines:",lines[part_index[k]])
        s = replace(lines[part_index[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        ss = ss[2]
        name_ = ss[ findfirst("name=", ss).stop + 1 : length(ss) ];
        PART[k].name = name_
        #println("PART[", k, "].name:", PART[k].name)

        #--- Node ---%
        index = 1;
        for i = part_index[k] : n
            if occursin.("*Node", lines[i]) #  strfind(lines{i}, '*Node') > 0
                index = i;
                break;
            end
        end

        nNode = 0;
        for i = index+1 : n
            if occursin.("*", lines[i]) #strfind(lines{i}, '*') > 0
                break;
            end
            nNode = nNode + 1;
        end
        PART[k].nNode = nNode;

        coordmat = zeros(nNode, 3);
        for i = 1 : nNode
          if occursin.("*", lines[index+i]) #strfind(lines{index+i}, '*') > 0
              break;
          end
          s = replace(lines[index+i], " " => "")  #s = replace(lines{index+i},' ', '');
          ss = split(s, ",", keepempty=false)  #split(s,',');
          coordmat[i,1] = parse(Float64, ss[2]) #str2double( ss[2]);
          coordmat[i,2] = parse(Float64, ss[3]) #str2double( ss[3]);
          coordmat[i,3] = parse(Float64, ss[4]) #str2double( ss[4]);
        end
        #PART[k].coordmat = coordmat
        PART[k].coordmat = coordmat'  # column major
  
        #--- Element ---%
        index = 1;
        for i = part_index[k] : n
          if occursin.("*Element", lines[i]) #strfind(lines{i}, '*Element') > 0
              index = i;
              break;
          end
        end
      
        nElement = 0;
        for i = index+1 : n
          if occursin.("*", lines[i]) #strfind(lines{i}, '*') > 0
              break;
          end
          nElement = nElement + 1;
        end
        PART[k].nElement = nElement;
        elementmat = zeros(Int, nElement, 8);
      
        for i = 1 : nElement
          if occursin.("*", lines[index+i]) #strfind(lines{index+i}, '*') > 0
              break;
          end
          s = replace(lines[index+i], " " => "")  #s = replace(lines{index+i},' ', '');
          ss = split(s, ",", keepempty=false)  #ss = split(s,',');
          for j = 1 : 8
              elementmat[i,j] = parse(Int, ss[1+j]) #str2num( ss{1+j});
          end
        end
        #PART[k].elementmat = elementmat
        PART[k].elementmat = elementmat' # column major

        #--- Nset ---%
        nset_index = [];
        for i = part_index[k] : n
          if occursin.("*Nset", lines[i]) && occursin.("generate", lines[i])  #contains(lines{i}, '*Nset') > 0  && contains(lines{i}, 'generate') > 0
              #nset_index = [nset_index i];
              push!(nset_index, i)
          end
          if occursin.("*End Part", lines[i])  #contains(lines{i}, '*End Part') > 0
              break
          end
        end
      
        nset_num = length(nset_index)
        #NSET_p(1:nset_num) = struct('name', "");
      
        for i = 1 : nset_num
            ns = NsetType("", "", 0, "", 0, zeros(Int,1))
            push!(PART[k].NSET, ns)
            
            s = replace(lines[nset_index[i]], " " => "")
            ss = split(s, ",", keepempty=false)
            ss = ss[2]
            name_ = ss[ findfirst("nset=", ss).stop + 1 : length(ss) ];
            PART[k].NSET[i].name = name_
            
            s = replace(lines[nset_index[i]+1], " " => "") #replace(lines{nset_index(i)+1},' ', '');
            ss = split( s,',');
            a = [j for j in parse(Int, ss[1]):parse(Int, ss[3]):parse(Int, ss[2])] #[parse(Int, ss[1]) : parse(Int, ss[3]) : parse(Int, ss[2])];
            PART[k].NSET[i].nodes = a;    
        end
        
        for i = part_index[k] : n
            if occursin.("*Solid Section", lines[i])  #strfind(lines{i}, '*Solid Section') > 0
              s = replace(lines[i], " " => "")
              ss = split(s, ",", keepempty=false)
              for j = 1 : length(ss)
                  sss = ss[j];
                  if occursin.("material=", sss)  #strfind(sss, 'material=') > 0
                    name_ = sss[ findfirst("material=", sss).stop + 1 : length(sss) ];
                    PART[k].material_name = name_;
                    break
                  end
              end
              break
            end
        end
  
    end


    #%--- Instance ---%
    instance_index = [];
    for i = 1 : n
      if occursin.("*Instance", lines[i])  #contains(lines{i}, '*Instance') > 0  
          #instance_index = [instance_index i];
          push!(instance_index, i)
      end
    end
    instance_num = length(instance_index);
    #INSTANCE(1:instance_num) = struct('name', "");
  
    INSTANCE = InstanceType[]

    for k = 1 : instance_num

        p = InstanceType("", "", 0, 0, Vector{String}[], 0, 0, 0, 0, zeros(Int,0), 
                         zeros(Int,0,0), zeros(Int,0,0), zeros(Int,0), zeros(Int,0,0), zeros(Int,0), zeros(Int,0) )
        push!(INSTANCE,p)

        s = replace(lines[instance_index[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        ss = ss[2]
        name_ = ss[ findfirst("name=", ss).stop + 1 : length(ss) ];
        INSTANCE[k].name = name_

        ss = split(s, ",", keepempty=false)
        ss = ss[3]
        name_ = ss[ findfirst("part=", ss).stop + 1 : length(ss) ];
        INSTANCE[k].part_name = name_;
  
        for i = 1 : nPart
            if PART[i].name == INSTANCE[k].part_name #strcmp(PART(i).name, INSTANCE(k).part_name )
              INSTANCE[k].part_id = i;
              break
            end
        end
  
        #INSTANCE(k).translate = {};
        j = 1;
        translate = []
        for i = instance_index[k]+1 : n
            if occursin.("*End Instance", lines[i]) #contains(lines{i}, '*End Instance') > 0  
                break
            end
  
            s = replace(lines[i], " " => "")
            push!(translate, s)
            j = j + 1;
        end
        INSTANCE[k].translate = translate
  
    end
  

    #%--- Nset ---%
    nset_index = [];
    for i = 1 : n
      if occursin.("*Nset", lines[i]) && occursin.("instance=", lines[i]) #contains(lines{i}, '*Nset') > 0  && contains(lines{i}, 'instance=') > 0
          #nset_index = [nset_index i];
          push!(nset_index, i)
      end
    end
  
    nset_num = length(nset_index);
    #NSET(1:nset_num) = struct('name', "");

    NSET = NsetType[]
  
    for k = 1 : nset_num

        p = NsetType("", "", 0, "", 0, zeros(Int,1))
        push!(NSET,p)

        s = replace(lines[nset_index[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        ss = ss[2]
        name_ = ss[ findfirst("nset=", ss).stop + 1 : length(ss) ];
        NSET[k].name = name_
       
        ss = split(s, ",", keepempty=false)
        ss = ss[3]
        name_ = ss[ findfirst("instance=", ss).stop + 1 : length(ss) ];
        NSET[k].instance_name = name_;
  
        for i = 1 : instance_num
            if NSET[k].instance_name == INSTANCE[i].name  #strcmp(NSET(k).instance_name, INSTANCE(i).name)
                  NSET[k].part_name = INSTANCE[i].part_name;
                  NSET[k].part_id = INSTANCE[i].part_id;
                  NSET[k].instance_id = i;
            end
        end
  
        a = [];
        index = nset_index[k];
  
        ss = split(s, ",", keepempty=false)
        if length(ss) == 4 && ss[4] == "generate"  #strcmp(ss{4},"generate")
  
             s = replace(lines[index+1], " " => "")
             ss = split(s, ",", keepempty=false)
             a = [j for j in parse(Int, ss[1]):parse(Int, ss[3]):parse(Int, ss[2])] #[str2num(ss{1}) : str2num(ss{3}) : str2num(ss{2})];
             
        else
  
            for i = index+1 : n
              if occursin.("*", lines[i])  #strfind(lines{i}, '*') > 0
                  break;
              end
      
              s = replace(lines[i], " " => "")
              ss = split(s, ",", keepempty=false)
              for j = 1 : length(ss)
                  if length(ss[j]) > 0
                      #a = [a  str2num( ss{j})];
                      push!(a, parse(Int, ss[j]))
                  end
              end
            end
  
        end
        NSET[k].nodes = a;
    end


    #%--- Elset ---%
    elset_index = [];
    for i = 1 : n
      if occursin.("*Elset", lines[i]) && occursin.("instance=", lines[i]) #contains(lines{i}, '*Nset') > 0  && contains(lines{i}, 'instance=') > 0
          push!(elset_index, i)
      end
    end
  
    elset_num = length(elset_index);
    ELSET = ELsetType[]
  
    for k = 1 : elset_num
         
        p = ELsetType("", "", 0, "", 0, zeros(Int,1))
        push!(ELSET,p)

        s = replace(lines[elset_index[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        ss = ss[2]
        name_ = ss[ findfirst("elset=", ss).stop + 1 : length(ss) ];
        ELSET[k].name = name_
       
        ss = split(s, ",", keepempty=false)
        ss = ss[3]
        if occursin.("instance=", ss)
            name_ = ss[ findfirst("instance=", ss).stop + 1 : length(ss) ];
            ELSET[k].instance_name = name_;    
        else
            ss = split(s, ",", keepempty=false)
            ss = ss[4]
            if occursin.("instance=", ss)
                name_ = ss[ findfirst("instance=", ss).stop + 1 : length(ss) ];
                ELSET[k].instance_name = name_;  
            end
        end
 
        for i = 1 : instance_num
            if ELSET[k].instance_name == INSTANCE[i].name  
                  ELSET[k].part_name = INSTANCE[i].part_name;
                  ELSET[k].part_id = INSTANCE[i].part_id;
                  ELSET[k].instance_id = i;
            end
        end

        a = [];
        index = elset_index[k];
  
        ss = split(s, ",", keepempty=false)
        if length(ss) == 4 && ss[4] == "generate"  
  
             s = replace(lines[index+1], " " => "")
             ss = split(s, ",", keepempty=false)
             a = [j for j in parse(Int, ss[1]):parse(Int, ss[3]):parse(Int, ss[2])] 

        elseif length(ss) == 5 && ss[3] == "internal" && ss[5] == "generate"
  
             s = replace(lines[index+1], " " => "")
             ss = split(s, ",", keepempty=false)
             a = [j for j in parse(Int, ss[1]):parse(Int, ss[3]):parse(Int, ss[2])] 
  
        elseif length(ss) == 4 && ss[3] == "internal"
  
            for i = index+1 : n
                if occursin.("*", lines[i])
                    break;
                end
        
                s = replace(lines[i], " " => "")
                ss = split(s, ",", keepempty=false)
                for j = 1 : length(ss)
                    if length(ss[j]) > 0
                        push!(a, parse(Int, ss[j]))
                    end
                end
            end
        
        end
        ELSET[k].elements = a;
  
    end
    #println("ELSET:", ELSET)

    #%--- Surface ---%
    surface_index = [];
    for i = 1 : n
      if occursin.("*Surface,", lines[i]) 
          push!(surface_index, i)
      end
    end
  
    surface_num = length(surface_index);
    SURFACE = SurfaceType[]

    for k = 1 : surface_num
         
        p = SurfaceType("", Vector{String}[], 0, zeros(Int,1))
        push!(SURFACE,p)

        s = replace(lines[surface_index[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        ss = ss[3]
        name_ = ss[ findfirst("name=", ss).stop + 1 : length(ss) ];
        SURFACE[k].name = name_
       
        index = surface_index[k];
        a = Int[];

        c = 1;
        for i = index+1 : n
            if occursin.("*", lines[i])
                break;
            end
    
            s = replace(lines[i], " " => "")
            ss = split(s, ",", keepempty=false)
            push!(SURFACE[k].elset_name, ss[1])
          
            for j = 1 : elset_num
                if SURFACE[k].elset_name[c] == ELSET[j].name
                    SURFACE[k].instance_id = ELSET[j].instance_id;
                    append!(a, ELSET[j].elements)
                end
            end
  
          c = c + 1;
        end
        SURFACE[k].elements = unique(sort(a));
           
    end  
    #println("SURFACE:", SURFACE)


    #%--- Global Model ---%
    nNode = 0;
    nElement = 0;
    coordmat = [];
    elementmat = [];
    for i = 1 : instance_num
        part_id = INSTANCE[i].part_id;
        coordmat_i = PART[part_id].coordmat;
        INSTANCE[i].node_offset = nNode;
        INSTANCE[i].element_offset = nElement;
        INSTANCE[i].nNode = PART[part_id].nNode;
        INSTANCE[i].nElement = PART[part_id].nElement;
        INSTANCE[i].elements = [k for k in 1 : INSTANCE[i].nElement]
  
        T =  Matrix(I,3,3)  #eye(3);
        for j = length(INSTANCE[i].translate): -1 : 1
            s = INSTANCE[i].translate[j];
            ss = split(s, ",", keepempty=false)
            if length(ss) == 3
                #offset_ = [parse(Float64, ss[1]) parse(Float64, ss[2]) parse(Float64, ss[3])];
                #coordmat_i = coordmat_i + ones(size(coordmat_i,1),1) * offset_;
                offset_ = [parse(Float64, ss[1]) parse(Float64, ss[2]) parse(Float64, ss[3])]'
                coordmat_i = coordmat_i + offset_ * ones(1,size(coordmat_i,2))
            elseif length(ss) == 7
                nv = [parse(Float64, ss[4]) - parse(Float64, ss[1])
                      parse(Float64, ss[5]) - parse(Float64, ss[2])
                      parse(Float64, ss[6]) - parse(Float64, ss[3]) ];
                nv = nv / norm(nv);
                n1 = nv[1];
                n2 = nv[2];
                n3 = nv[3];
                d = parse(Float64, ss[7]) / 180.0 * pi;
                T = [n1*n1*(1-cos(d))+cos(d)  n1*n2*(1-cos(d))-n3*sin(d)  n1*n3*(1-cos(d))+n2*sin(d)
                     n1*n2*(1-cos(d))+n3*sin(d)  n2*n2*(1-cos(d))+cos(d)  n2*n3*(1-cos(d))-n1*sin(d)
                     n1*n3*(1-cos(d))-n2*sin(d)  n2*n3*(1-cos(d))+n1*sin(d)  n3*n3*(1-cos(d))+cos(d) ];
                #coordmat_i = (T * coordmat_i')';
                coordmat_i = T * coordmat_i;
            end
        end

        elementmat_i = PART[part_id].elementmat .+ nNode;
  
        if i == 1
            coordmat = coordmat_i;
            elementmat = elementmat_i;
        else
            #coordmat = [coordmat; coordmat_i];
            #elementmat = [elementmat; elementmat_i];
            coordmat = [coordmat  coordmat_i];
            elementmat = [elementmat  elementmat_i];
        end
        
        nNode = nNode + PART[part_id].nNode;
        nElement = nElement .+ PART[part_id].nElement;
    end
  

    #%--- Amplitude ---%
    amplitude_index = [];
    for i = 1 : n
      if occursin.("*Amplitude", lines[i]) #contains(lines{i}, '*Amplitude') > 0  
          push!(amplitude_index, i)
      end
    end
  
    amplitude_num = length(amplitude_index);
    #AMPLITUDE(1:amplitude_num) = struct('name', "", 'time', [], 'value', []);

    AMPLITUDE = AmplitudeType[]
   
    for k = 1 : amplitude_num

        p = AmplitudeType("", zeros(Float64,1), zeros(Float64,1))
        push!(AMPLITUDE,p)        

        s = replace(lines[amplitude_index[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        ss = ss[2]
        name_ = ss[ findfirst("name=", ss).stop + 1 : length(ss) ];
        AMPLITUDE[k].name = name_

        index = amplitude_index[k];  
        for i = index+1 : n
          if occursin.("*", lines[i]) #strfind(lines{i}, '*') > 0
              break;
          end
          
          s = replace(lines[i], " " => "")
          ss = split(s, ",", keepempty=false)
          time_ = []
          value_ = []
          for j = 1 : Int(length(ss)/2)
              push!(time_, parse(Float64, ss[j*2-1]))
              push!(value_, parse(Float64, ss[j*2]))
              #println("time_:", parse(Float64, ss[j*2-1]))
          end
          AMPLITUDE[k].time = time_
          AMPLITUDE[k].value = value_
  
        end
  
    end

    
    #%--- Material ---%
    material_index = [];
    for i = 1 : n
      if occursin.("*Material", lines[i]) #contains(lines{i}, '*Material') > 0  
          push!(material_index, i)
      end
    end
  
    material_num = length(material_index);
    #METERIAL(1:material_num) = struct('name', "");

    MATERIAL = MaterialType[]
  
    for k = 1 : material_num

        p = MaterialType("", 0.0, 0.0, 0.0, zeros(Float64, 0,0), zeros(Float64, 0), 0, 0.0, zeros(Float64, 0,0), 0.0, zeros(1,1))
        push!(MATERIAL,p)        

        s = replace(lines[material_index[k]], " " => "")
        ss = split(s, ",", keepempty=false)
        ss = ss[2]
        name_ = ss[ findfirst("name=", ss).stop + 1 : length(ss) ];
        MATERIAL[k].name = name_
        #MATERIAL[k].fracture_flag = 0;
        #MATERIAL[k].failure_stress = [];
  
        index = material_index[k];
        plastic_index = 1;
        ductile_index = 1;
  
        for i = index+1 : n
          if occursin.("*Material", lines[i]) #strfind(lines{i}, '*Material') > 0
              break;
          end
          if occursin.("**", lines[i]) #strfind(lines{i}, '**') > 0
              break;
          end
  
          if occursin.("*Density", lines[i]) #strfind(lines{i}, '*Density') > 0
              s = replace(lines[i+1], " " => "")
              ss = split(s, ",", keepempty=false)
              MATERIAL[k].density = parse(Float64, ss[1]) #str2num( ss{1} );
          end
  
          if occursin.("*Elastic", lines[i]) #strfind(lines{i}, '*Elastic') > 0
              s = replace(lines[i+1], " " => "")
              ss = split(s, ",", keepempty=false)
              MATERIAL[k].young = parse(Float64, ss[1]) #str2num( ss{1} );
              MATERIAL[k].poisson = parse(Float64, ss[2]) #str2num( ss{2} );
          end
  
          if occursin.("*Plastic", lines[i]) #strfind(lines{i}, '*Plastic') > 0
              plastic_index = i;
          end
  
          if occursin.("*Damage Initiation", lines[i]) && occursin.("criterion=DUCTILE", lines[i]) #contains(lines{i}, '*Damage Initiation')  && contains(lines{i}, 'criterion=DUCTILE')
              ductile_index = i
              MATERIAL[k].fracture_flag = 1;
          end
          
          if occursin.("*Tensile Failure", lines[i]) #strfind(lines{i}, '*Tensile Failure') > 0
              s = replace(lines[i+1], " " => "")
              ss = split(s, ",", keepempty=false)
              MATERIAL[k].failure_stress = parse(Float64, ss[1]) #str2num( ss{1} );
              MATERIAL[k].fracture_flag = 1;
          end
  
        end
  
        plastic = [];
        if plastic_index > material_index[k]
            for i = plastic_index+1 : n
              if occursin.("*", lines[i]) #strfind(lines{i}, '*') > 0
                  break;
              end
              s = replace(lines[i], " " => "")
              ss = split(s, ",", keepempty=false)
              #println("ss:", ss)
              #plastic = [plastic; str2num(ss{1})  str2num(ss{2})];
              push!(plastic, [parse(Float64, ss[1])  parse(Float64, ss[2])])
            end
            #MATERIAL[k].plastic = plastic;

            m_ = zeros(Float64, length(plastic), 2)
            for i = 1 : length(plastic)
                m_[i,1] = plastic[i][1]
                m_[i,2] = plastic[i][2]
            end
            MATERIAL[k].plastic = m_;
        end
        #println("plastic:", MATERIAL[k].plastic)

        Hd = Float64[]
        for i = 1 : size(MATERIAL[k].plastic,1) - 1
            v = (MATERIAL[k].plastic[i+1,1] - MATERIAL[k].plastic[i,1]) / (MATERIAL[k].plastic[i+1,2] - MATERIAL[k].plastic[i,2])
            append!(Hd, v)
        end
        MATERIAL[k].Hd = Hd
        #println("Hd:", MATERIAL[k].Hd)
        
        ductile = [];
        if ductile_index > material_index[k]
            for i = ductile_index+1 : n
              if occursin.("*", lines[i]) #strfind(lines{i}, '*') > 0
                  break;
              end
              s = replace(lines[i], " " => "")
              ss = split(s, ",", keepempty=false)
              #ductile = [ductile; str2num(ss{1})  str2num(ss{2})  str2num(ss{3})];
              push!(ductile, [parse(Float64, ss[1])  parse(Float64, ss[2])  parse(Float64, ss[3])])
            end
            #MATERIAL[k].ductile = ductile;

            m_ = zeros(Float64, length(ductile), 3)
            for i = 1 : length(ductile)
                m_[i,1] = ductile[i][1]
                m_[i,2] = ductile[i][2]
                m_[i,3] = ductile[i][3]
            end
            MATERIAL[k].ductile = m_;
        end
               
    end

    
    element_material = Int[];
    element_instance = Int[];
    for i = 1 : length(INSTANCE) 
        part_id = INSTANCE[i].part_id;
         for j = 1 : length(MATERIAL)
             if PART[part_id].material_name == MATERIAL[j].name #strcmp(PART(part_id).material_name, MATERIAL(j).name)
                  PART[part_id].material_id = j;
                  INSTANCE[i].material_id = j;
             end
         end
         mat_id = ones(Int, PART[part_id].nElement, 1) .* PART[part_id].material_id;
         #element_material = [element_material; mat_id];
         append!(element_material, mat_id)
  
         instance_id = ones(Int, PART[part_id].nElement, 1) .* i;
         #element_instance = [element_instance; instance_id];
         append!(element_instance, instance_id)
    end
  
    
    #%--- Step ---%
    d_time = 0.;
    end_time = 0.;
    for i = 1 : n
      if occursin.("*Dynamic, Explicit", lines[i])  #contains(lines{i}, '*Dynamic, Explicit') > 0 
          s = replace(lines[i+1], " " => "")
          ss = split(s, ",", keepempty=false) 
          d_time = parse(Float64, ss[1]) #str2double(ss{1});
          end_time = parse(Float64, ss[2]) #str2double(ss{2});
          break
      end
    end
  
    #%--- Mass scaling ---%
    mass_scaling = 1.;
    for i = 1 : n
      if occursin.("*Fixed Mass Scaling", lines[i])  #contains(lines{i}, '*Fixed Mass Scaling') > 0  
          s = replace(lines[i], " " => "")
          ss = split(s, ",", keepempty=false) 
          ss = ss[2];
          v_ = ss[ findfirst("factor=", ss).stop + 1 : length(ss) ];
          mass_scaling = parse(Float64, v_) #str2double(v_);
          break
      end
    end

    
    #%--- BC ---%
    bc_index = [];
    for i = 1 : n
      if occursin.("*Boundary", lines[i]) #contains(lines{i}, '*Boundary') > 0  
          push!(bc_index, i) 
      end
    end
  
    bc_num = length(bc_index);
    BC = BCType[]
  
    for k = 1 : bc_num

        p = BCType("", [], [], "", AmplitudeType("", zeros(Float64,1), zeros(Float64,1)))
        push!(BC,p)        

        index = bc_index[k];
        
        s = replace(lines[index], " " => "")
        ss = split(s, ",", keepempty=false) 
        if length(ss) == 2 && occursin.("amplitude=", ss[2]) #contains(ss{2}, 'amplitude=')
            ss = ss[2];
            name_ = ss[ findfirst("amplitude=", ss).stop + 1 : length(ss) ];
            BC[k].amp_name = name_;
            for j = 1 : amplitude_num
                  if AMPLITUDE[j].name == name_ #strcmp(AMPLITUDE(j).name, name_)
                      BC[k].amplitude = AMPLITUDE[j];
                      break
                  end
            end
        end
  
        
        for i = index+1 : n
          if occursin.("*Boundary", lines[i]) #strfind(lines{i}, '*Boundary') > 0
              break;
          end
          if occursin.("**", lines[i]) #strfind(lines{i}, '**') > 0
              break;
          end
  
          s = replace(lines[i], " " => "")
          ss = split(s, ",", keepempty=false) 
          BC[k].Nset_name = ss[1];
  
          nodes = [];
          if occursin.(".", BC[k].Nset_name) #strfind(BC[k].Nset_name, '.') > 0
              sss = split(BC[k].Nset_name,'.', keepempty=false) 
              instance_name = sss[1];
              nset_name = sss[2];
              instance_id = 0;
              part_id = 0;
              for j = 1 : instance_num
                  if INSTANCE[j].name == instance_name #strcmp(INSTANCE[j].name, instance_name)
                      instance_id = j;
                      part_id = INSTANCE[j].part_id;
                      break
                  end
              end
  
              for j = 1 : length(PART[part_id].NSET)
                  if PART[part_id].NSET[j].name == nset_name #strcmp(PART[part_id].NSET[j].name, nset_name)
                      nodes_j = PART[part_id].NSET[j].nodes .+ INSTANCE[instance_id].node_offset;
                      #nodes = [nodes  nodes_j];
                      append!(nodes, nodes_j)
                      break
                  end
              end
  
          else 
              for j = 1 : nset_num
                  if BC[k].Nset_name == NSET[j].name #strcmp(BC[k].Nset_name, NSET[j].name)
                      nodes_j = NSET[j].nodes .+ INSTANCE[NSET[j].instance_id].node_offset;
                      #nodes = [nodes nodes_j];
                      append!(nodes, nodes_j)
                  end
              end
          end
  
          dof = Int[];
          if length(ss) == 2 && occursin.("ENCASTRE", ss[2]) #strcmp(ss{2}, 'ENCASTRE')
              #dof = [ nodes*3 .- 2, nodes*3 .- 1,  nodes*3  ];
              append!(dof, nodes*3 .- 2)
              append!(dof, nodes*3 .- 1)
              append!(dof, nodes*3 )
              #BC[k].dof = sort!(dof);
              push!(BC[k].dof, dof)
              BC[k].value = [0.];
          elseif length(ss) == 3
              ii = parse(Int, ss[2]) #str2num(ss[2]);
              dir = parse(Int,ss[3]);
              if dir <= 3
                  dof = nodes*3 .- (3-dir);
#                  BC[k].dof{ii} = dof;
#                  BC[k].value{ii} = 0;
                  push!(BC[k].dof, dof)
                  push!(BC[k].value, 0.)
              end
          elseif length(ss) == 4
              ii = parse(Int,ss[2]);
              dir = parse(Int,ss[3]);
              value = parse(Float64,ss[4]);
              if dir <= 3
                  dof = nodes*3 .- (3-dir);
 #                 BC[k].dof{ii} = dof;
 #                 BC[k].value{ii} = value;
                   push!(BC[k].dof, dof)
                   push!(BC[k].value, value)
              end
  
          end
  
        end
    
    end
  

    #%--- Initial Conditions ---%
    ic_index = [];
    for i = 1 : n
      if occursin.("*Initial Conditions", lines[i]) #contains(lines{i}, '*Initial Conditions') > 0  
          push!(ic_index, i)
      end
    end
  
    ic_num = length(ic_index);
    IC = ICType[]

    for k = 1 : ic_num

        p = ICType("", "", [], [])
        push!(IC,p)        

        index = ic_index[k];

        s = replace(lines[index], " " => "")
        ss = split(s, ",", keepempty=false) 
        ss = ss[2];
        type_ = ss[ findfirst("type=", ss).stop + 1 : length(ss) ];
        IC[k].type = type_;
        
        for i = index+1 : n
          if occursin.("*Initial Conditions", lines[i]) # strfind(lines{i}, '*Initial Conditions') > 0
              break;
          end
          if occursin.("**", lines[i]) #strfind(lines{i}, '**') > 0
              break;
          end
  
          s = replace(lines[i], " " => "")
          ss = split(s, ",", keepempty=false) 
          IC[k].Nset_name = ss[1];
  
          nodes = Int[];
          if occursin.(".", IC[k].Nset_name)  #strfind(IC[k].Nset_name, '.') > 0
              sss = split(IC[k].Nset_name,'.', keepempty=false) 
              instance_name = sss[1];
              nset_name = sss[2];
              instance_id = 0;
              part_id = 0;
              for j = 1 : instance_num
                  if INSTANCE[j].name == instance_name #strcmp(INSTANCE[j].name, instance_name)
                      instance_id = j;
                      part_id = INSTANCE[j].part_id;
                      break
                  end
              end
  
              for j = 1 : length(PART[part_id].NSET)
                  if PART[part_id].NSET[j].name == nset_name #strcmp(PART[part_id].NSET[j].name, nset_name)
                      nodes_j = PART[part_id].NSET[j].nodes .+ INSTANCE[instance_id].node_offset;
                      nodes = nodes_j;
                      break
                  end
              end
  
          else 
              for j = 1 : nset_num
                  if IC[k].Nset_name == NSET[j].name #strcmp(IC[k].Nset_name, NSET[j].name)
                      nodes_j = NSET[j].nodes .+ INSTANCE[NSET[j].instance_id].node_offset;
                      nodes = nodes_j;
                      break
                  end
              end
          end
  
          #dof = [];
          dir = parse(Int,ss[2]) #str2num(ss{2});
          dof = nodes*3 .- (3-dir);
          #println("dof:", dof)
  
          v = parse(Float64,ss[3])  # str2num(ss{3});
          
          #IC[k].dof = dof;
          #IC[k].value = v;
          push!(IC[k].dof, dof)
          push!(IC[k].value, v)
  
        end
  
    end
  

    #%--- Interactions (Contact) ---%
    contact_flag = 0;
    for i = 1 : n
      if occursin.("*Contact", lines[i]) # contains(lines{i}, '*Contact') > 0   %'*Contact Inclusions'
          contact_flag = 1;
          break
      end
    end

    for i = 1 : n
        if occursin.("*Contact Inclusions", lines[i]) && occursin.("HAKAIoption=self-contact", lines[i])
            contact_flag = 2;
            break
        end
      end
  

    cp_index = [];
    for i = 1 : n
      if occursin.("*Contact Pair,", lines[i]) #contains(lines{i}, '*Initial Conditions') > 0  
          push!(cp_index, i)
      end
    end
  
    cp_num = length(cp_index);
    CP = CPType[]
    
    for k = 1 : cp_num

        p = CPType("", "", "", 0, 0, Int[], Int[],
                   zeros(Int,0,0), zeros(Int,0,0), zeros(Int,0), zeros(Int,0), zeros(Int,0), zeros(Int,0))
        push!(CP,p)     

        index = cp_index[k];
        s = replace(lines[index], " " => "")
        ss = split(s, ",", keepempty=false)
        ss = ss[4]
        name_ = ss[ findfirst("cpset=", ss).stop + 1 : length(ss) ];
        CP[k].name = name_

        s = replace(lines[index+1], " " => "")
        ss = split(s, ",", keepempty=false)
        CP[k].surface_name_1 = ss[1];
        CP[k].surface_name_2 = ss[2];
          
        for j = 1 : surface_num
            if CP[k].surface_name_1 == SURFACE[j].name
                CP[k].instance_id_1 = SURFACE[j].instance_id;
                CP[k].elements_1 = SURFACE[j].elements;
            end
            if CP[k].surface_name_2 == SURFACE[j].name
                CP[k].instance_id_2 = SURFACE[j].instance_id;
                CP[k].elements_2 = SURFACE[j].elements;
            end
        end
  
    end
    #println("CP:", CP)


    MODEL = ModelType(PART, INSTANCE, NSET, ELSET, SURFACE, AMPLITUDE, MATERIAL, BC, IC, CP, 
                      nNode, coordmat, nElement, elementmat, element_material, element_instance,
                      d_time, end_time, mass_scaling, contact_flag)
    

    return MODEL  #PART

end