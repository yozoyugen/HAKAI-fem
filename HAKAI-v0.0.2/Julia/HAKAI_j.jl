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
using StaticArrays
using Printf
using Plots
using Random
using Quadmath
using Base.Threads
using FLoops
using CUDA

include("./readInpFile_j.jl")


mutable struct IntegDataType
    index::Int
    integ_stress::Array{Float64,2}
    integ_strain::Array{Float64,2}
    integ_plastic_strain::Array{Float64,2}
    integ_eq_plastic_strain::Array{Float64,1}
    integ_triax_stress::Array{Float64,1}
end

mutable struct NodeDataType
    node_stress::Array{Float64,2}
    node_strain::Array{Float64,2}
    node_plastic_strain::Array{Float64,2}
    node_eq_plastic_strain::Array{Float64,1}
    node_mises_stress::Array{Float64,1}
    node_triax_stress::Array{Float64,1}
end

mutable struct ContactValue
    A::Array{Float64,2}
    b::Array{Float64,1}
    n::Array{Float64,1}
    x::Array{Float64,1}

    p::Array{Float64,1}
    q0::Array{Float64,1}
    q1::Array{Float64,1}
    q2::Array{Float64,1}
    v1::Array{Float64,1}
    v2::Array{Float64,1}
    
    ve::Array{Float64,1}
    vs::Array{Float64,1}
    f::Array{Float64,1}
    fric::Array{Float64,1}
    v::Array{Float64,1}
end

mutable struct ContactTriangle
    c_nodes_i::Array{Int,1}
    c_nodes_j::Array{Int,1}
    c_triangles_eleid::Array{Int,1}
    c_triangles::Array{Int,2}
    young::Float64
end


function hakai(fname)

    Nth = Threads.nthreads()
    println("nthreads:", Nth)

    rname = @sprintf("julia_bug_report.txt")
    bug_report = open(rname,"w")

    MODEL = readInpFile( fname )
    #println("MODEL:", MODEL)
    #println("MODEL.PART:", MODEL.PART)
    #println("MODEL.PART[1]:", MODEL.PART[1])
    #println("MODEL.PART[1].name:", MODEL.PART[1].name)

    nNode = MODEL.nNode;
    println("nNode:", nNode)
    coordmat = MODEL.coordmat # column major
    nElement = MODEL.nElement;
    println("nElement:", nElement)
    elementmat = MODEL.elementmat # column major
    element_material = MODEL.element_material;
    element_instance = MODEL.element_instance;
    contact_flag = MODEL.contact_flag;
    println("contact_flag:", contact_flag)

    gpu_contact_flag = 0
    println("gpu_contact_flag:", gpu_contact_flag)

    
    flag_fracture = 0;
    
    mass_scaling = MODEL.mass_scaling;
    println("mass_scaling:", mass_scaling)
    d_time = MODEL.d_time * sqrt(mass_scaling);
    end_time = MODEL.end_time;
    time_num = end_time / d_time;
    println("time_num:", time_num)

    div_d_time = 1.0 / d_time
    div_d_time2 = 1.0 / d_time^2

    force_dof = Int[];
    force_v = Float64[];

#=    c_pair = Vector{Vector{Int}}();
    if contact_flag == 1
        ni = length(MODEL.INSTANCE);
        #%ni = 2;
        for i = 1 : ni
            for j = 1 : ni
                if i != j
                    #c_pair = [c_pair; i j];
                    push!(c_pair, [i,j])
                end
            end
        end    
        #c_pair
        println("c_pair:", c_pair)
    end =#


    #%--- Material property 
    for i = 1 : length(MODEL.MATERIAL)
        young = MODEL.MATERIAL[i].young;
        poisson = MODEL.MATERIAL[i].poisson;
        G = young / 2. / (1.0+poisson);
        MODEL.MATERIAL[i].G = G;
    
        #%Dmat = zeros(6,6);
        d1 = (1.0-poisson);
        d2 = poisson;
        d3 = (1.0-2.0*poisson)/2.0;
        Dmat_i = young / (1.0+poisson) / (1.0-2.0*poisson) *
                          [ d1 d2 d2 0  0  0
                            d2 d1 d2 0  0  0
                            d2 d2 d1 0  0  0
                            0  0  0  d3 0  0
                            0  0  0  0  d3 0
                            0  0  0  0  0  d3] ; #@SMatrix
        MODEL.MATERIAL[i].Dmat = Dmat_i;

        if length( MODEL.MATERIAL[i].failure_stress ) > 0
            flag_fracture = 1;
            #MODEL.MATERIAL[i].failure_stress
        end

        if size( MODEL.MATERIAL[i].ductile, 1) > 0
            flag_fracture = 1;        
        end
        
        #println("MATERIAL[", i, "]:", MODEL.MATERIAL[i] )
    end
   

    # set variable
    fn = nNode * 3;
    integ_num = 8;

    # shape func
    Pusai_mat = cal_Pusai_hexa(integ_num);
    #println("Pusai_mat:", Pusai_mat)

    elementVolume = zeros(nElement);
    for e = 1 : nElement
        #e_position = coordmat[elementmat[e,:], :];
        #e_position = coordmat[elementmat[:,e], :];
        e_position = coordmat[:, elementmat[:,e]];
        V = 0.;
        for i = 1 : integ_num
            #J = Pusai_mat[:,:,i] * e_position;
            #J = Pusai_mat[:,:,i] * e_position';
            J = Pusai_mat[i] * e_position';
            #V = V + det(J);
            V = V + my3det(J);
        end
        elementVolume[e] = V;
    end
    #elementVolume.';
    #print(bug_report, @sprintf("time=%.16e, V=%.16e\n", 0.0, sum(elementVolume)) ) 


    diag_M = zeros(fn);
    diag_C = zeros(fn);
    for i = 1 : nElement
        density = MODEL.MATERIAL[ element_material[i] ].density;
        node_mass = density * elementVolume[i] / 8.0;
        #diag_M[ (elementmat[i,:] .- 1)*3 .+ 1 ] .+= node_mass;
        #diag_M[ (elementmat[i,:] .- 1)*3 .+ 2 ] .+= node_mass;
        #diag_M[ (elementmat[i,:] .- 1)*3 .+ 3 ] .+= node_mass;
        diag_M[ (elementmat[:,i] .- 1)*3 .+ 1 ] .+= node_mass;
        diag_M[ (elementmat[:,i] .- 1)*3 .+ 2 ] .+= node_mass;
        diag_M[ (elementmat[:,i] .- 1)*3 .+ 3 ] .+= node_mass;
    end

    diag_M .= diag_M * mass_scaling;

    C = 0.0;
    diag_C .= diag_M * C;

    #for i = 1 : fn
    #    print(bug_report, @sprintf("%d, %.16e\n", i, diag_M[i] ) );
    #end


    position = copy(coordmat);
    disp = zeros(fn);
    disp_new = zeros(fn);
    disp_pre = zeros(fn);
    d_disp = zeros(fn);
    velo = zeros(fn);

    #%--- Initial Velocity ---%
    for i = 1 : length(MODEL.IC)
        #println("IC:", MODEL.IC[i])
        for j = 1 : length(MODEL.IC[i].dof)
            disp_pre[ MODEL.IC[i].dof[j] ] .= -MODEL.IC[i].value[j] * d_time
            velo[ MODEL.IC[i].dof[j] ] .= MODEL.IC[i].value[j]
        end
    end
    #println("disp_pre:", disp_pre)
    #println("velo:", velo)

    #%--- Contact ---%
    instance_pair = Vector{Vector{Int}}()
    cp_index = Int[];
    CV = ContactValue(zeros(3,3), zeros(3), zeros(3), zeros(3), 
                      zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3),
                      zeros(3), zeros(3), zeros(3), zeros(3), zeros(3) )
    all_exterior_flag = 0  # ->  0: CP based local contact,   1: all surface contact

    if contact_flag >= 1
        
        for i = 1 : length(MODEL.INSTANCE)  
    
            faces, faces_eleid, sorted_faces = get_element_face(MODEL, i);
            #println("faces:", faces)
            #println("faces_eleid:", faces_eleid)
            #println("sorted_faces:", sorted_faces)
            MODEL.INSTANCE[i].surfaces = faces;
            MODEL.INSTANCE[i].surfaces_eleid = faces_eleid;
            MODEL.INSTANCE[i].sorted_surfaces = sorted_faces;
    
            #=
            array_element =  [j for j in 1:MODEL.INSTANCE[i].nElement] #[1 : MODEL.INSTANCE[i].nElement];
            c_triangles, c_triangles_eleid, c_nodes = get_surface_triangle(MODEL, i, array_element);
            MODEL.INSTANCE[i].c_triangles = c_triangles;
            MODEL.INSTANCE[i].c_triangles_eleid = c_triangles_eleid;
            MODEL.INSTANCE[i].c_nodes = c_nodes;  =#
    
        end

        if length(MODEL.CP) == 0  #% ALL EXTERIOR

            all_exterior_flag = 1
            ni = length(MODEL.INSTANCE);
            CP = CPType[]

            if ni > 1
                c = 1;
                for i = 1 : ni
                    js = i+1;
                    #js = i;  #% -> include self-contact
                    if contact_flag == 2
                        js = i;
                    end
                    
                    for j = js : ni  
                        #% c_pair = [c_pair; i j];
                        #% c_pair = [c_pair; j i];

                        p = CPType("", "", "", 0, 0, Int[], Int[],
                                   zeros(Int,0,0), zeros(Int,0,0), zeros(Int,0), zeros(Int,0), zeros(Int,0), zeros(Int,0))
                        push!(CP,p)
                        CP[c].instance_id_1 = i;
                        CP[c].instance_id_2 = j;
                        CP[c].elements_1 = [k for k in 1 : MODEL.INSTANCE[i].nElement]
                        CP[c].elements_2 = [k for k in 1 : MODEL.INSTANCE[j].nElement]
                        #println("CP[", c, "]:", CP[c])
                        println("CP[", c, "]:", i, ", ", j)
                        c = c+1;
                    end
                end    

            else #% ni == 1
                p = CPType("", "", "", 0, 0, Int[], Int[],
                           zeros(Int,0,0), zeros(Int,0,0), zeros(Int,0), zeros(Int,0), zeros(Int,0), zeros(Int,0))
                push!(CP,p)
                CP[1].instance_id_1 = 1;
                CP[1].instance_id_2 = 1;
                CP[1].elements_1 = [k for k in 1 : MODEL.INSTANCE[1].nElement]
                CP[1].elements_2 = [k for k in 1 : MODEL.INSTANCE[1].nElement]
            end

            MODEL.CP = CP

        else

        end


        for i = 1 : length(MODEL.CP)
            array_element =  [j for j in 1:MODEL.INSTANCE[ MODEL.CP[i].instance_id_1 ].nElement] 
            c_triangles, c_triangles_eleid, c_nodes = get_surface_triangle(MODEL.INSTANCE[ MODEL.CP[i].instance_id_1 ], array_element,  MODEL.CP[i].elements_1);
            
            MODEL.CP[i].c_triangles_1 = c_triangles;
            MODEL.CP[i].c_triangles_eleid_1 = c_triangles_eleid;
            MODEL.CP[i].c_nodes_1 = c_nodes;

            array_element =  [j for j in 1:MODEL.INSTANCE[ MODEL.CP[i].instance_id_2 ].nElement] 
            c_triangles, c_triangles_eleid, c_nodes = get_surface_triangle(MODEL.INSTANCE[ MODEL.CP[i].instance_id_2 ], array_element,  MODEL.CP[i].elements_2);
            MODEL.CP[i].c_triangles_2 = c_triangles;
            MODEL.CP[i].c_triangles_eleid_2 = c_triangles_eleid;
            MODEL.CP[i].c_nodes_2 = c_nodes;

            #println("MODEL.CP[i]:", MODEL.CP[i])
        end

    
        for cc = 1 : length(MODEL.CP)
            if MODEL.CP[cc].instance_id_1 == MODEL.CP[cc].instance_id_2
                    #instance_pair = [instance_pair; CP(cc).instance_id_i  CP(cc).instance_id_j];
                    #cp_index = [cp_index cc];
                push!(instance_pair, [MODEL.CP[cc].instance_id_1, MODEL.CP[cc].instance_id_2])
                append!(cp_index, cc)
            else
                    #instance_pair = [instance_pair; CP(cc).instance_id_i  CP(cc).instance_id_j];
                    #instance_pair = [instance_pair; CP(cc).instance_id_j  CP(cc).instance_id_i];
                    #cp_index = [cp_index cc cc];
                push!(instance_pair, [MODEL.CP[cc].instance_id_1, MODEL.CP[cc].instance_id_2])
                push!(instance_pair, [MODEL.CP[cc].instance_id_2, MODEL.CP[cc].instance_id_1])
                append!(cp_index, cc)
                append!(cp_index, cc)
            end
        end         
        println("instance_pair:", instance_pair)

        CT = ContactTriangle[]
        # CP -> [1,2], ...
        # instance_pair, CT -> [1,2], [2,1], ...

        for c = 1 : length(cp_index)

            p = ContactTriangle(zeros(Int,0), zeros(Int,0), zeros(Int,0), zeros(Int,0,0), 0.0)
            push!(CT,p)

            cc = cp_index[c]
    
            #% j -> triangle,   i -> point
            i_instance = instance_pair[c][1]
            j_instance = instance_pair[c][2]
            #println("i_instance:", i_instance, " j_instance:", j_instance)
                
            young = MODEL.MATERIAL[ MODEL.INSTANCE[j_instance].material_id ].young;
            c_nodes_i = zeros(Int,0)  # -> 1 alloc
            c_nodes_j = zeros(Int,0)  # -> 1 alloc
            c_triangles_eleid = zeros(Int,0)  # -> 1 alloc
            c_triangles = zeros(Int,0,0)  # -> 1 alloc
    
            if MODEL.CP[cc].instance_id_1 == i_instance   # -> 4~5 alloc
                c_nodes_i = MODEL.CP[cc].c_nodes_1 .+ MODEL.INSTANCE[i_instance].node_offset
                c_nodes_j = MODEL.CP[cc].c_nodes_2 .+ MODEL.INSTANCE[j_instance].node_offset
                c_triangles = MODEL.CP[cc].c_triangles_2 .+ MODEL.INSTANCE[j_instance].node_offset
                c_triangles_eleid = MODEL.CP[cc].c_triangles_eleid_2 .+ MODEL.INSTANCE[j_instance].element_offset
            else
                c_nodes_i = MODEL.CP[cc].c_nodes_2 .+ MODEL.INSTANCE[i_instance].node_offset
                c_nodes_j = MODEL.CP[cc].c_nodes_1 .+ MODEL.INSTANCE[j_instance].node_offset
                c_triangles = MODEL.CP[cc].c_triangles_1 .+ MODEL.INSTANCE[j_instance].node_offset
                c_triangles_eleid = MODEL.CP[cc].c_triangles_eleid_1 .+ MODEL.INSTANCE[j_instance].element_offset
            end

            CT[c].c_nodes_i = c_nodes_i
            CT[c].c_nodes_j = c_nodes_j
            CT[c].c_triangles = c_triangles
            CT[c].c_triangles_eleid = c_triangles_eleid
            CT[c].young = young

            #println("CT[", c, "]:", CT[c])  # -> Match with Matlab
        end
  
    
        println("Contact set")
    end 


    elementSize = zeros(nElement,3);
    for e = 1 : nElement
        #e_pos = coordmat[elementmat[e,:], :];
        #L1 = my3norm( e_pos[1,:] - e_pos[2,:] );
        #L2 = my3norm( e_pos[1,:] - e_pos[4,:] );
        #L3 = my3norm( e_pos[1,:] - e_pos[5,:] );
        #e_pos = coordmat[elementmat[:,e], :];
        e_pos = coordmat[:,elementmat[:,e]];
        L1 = my3norm( e_pos[:,1] - e_pos[:,2] );
        L2 = my3norm( e_pos[:,1] - e_pos[:,4] );
        L3 = my3norm( e_pos[:,1] - e_pos[:,5] );
        elementSize[e,:] = [L1, L2, L3];
    end
    elementMinSize = minimum(elementSize)
    println("elementMinSize:", elementMinSize)
    elementMaxSize = maximum(elementSize)
    println("elementMaxSize:", elementMaxSize)

    d_max = 0.0
    #d_node = zeros(nNode)  #@MVector -> huge increase alloc
    #d_node_pre = zeros(nNode)
    d_disp_norm = zeros(nNode)
    

    #%--- Variable ---%
    external_force = zeros(fn);

    #c_force3 = zeros(fn);
    #c_force3 = zeros(fn, Nth);
    #c_force3 = zeros(Float32, fn, Nth);
    c_force3 = zeros(Float128, fn, Nth);

    Q = zeros(fn);
    Qe = zeros(24, nElement);
    #Q = [Threads.Atomic{Float64}(0) for _ in 1:fn]
            #Barray = zeros(6,24,integ_num) # -> even
        #Barray = [ @MMatrix zeros(6,24) for i=1:8]
        #BVarray = [ @MMatrix zeros(6,24) for i=1:8]
        #detJarray = @MVector zeros(integ_num)
        #Dmat = @MMatrix zeros(6,6)
        #e_position = @MMatrix zeros(3,8)

    #integ_stress = zeros( nElement * integ_num, 6);  # row major
    integ_stress = zeros( 6, nElement * integ_num); # column major
    integ_strain = zeros( 6, nElement * integ_num);
    integ_plastic_strain = zeros( 6, nElement * integ_num);
    integ_eq_plastic_strain = zeros( nElement * integ_num);
    integ_triax_stress = zeros( nElement * integ_num);
    element_flag = ones(Int, nElement);
    integ_flag = ones(Int, nElement, integ_num);

    integ_yield_stress = zeros( nElement * integ_num);
    for i = 1 : nElement
        pp = MODEL.MATERIAL[element_material[i]].plastic;
        #println("pp:", pp)
        #println("pp[1,1]:", pp[1,1])
        #println("pp[1][1]:", pp[1][1])
        if length( pp ) >0
            integ_yield_stress[(1:integ_num) .+ (i-1)*integ_num ] .= pp[1,1];
        end
    end
    #println("integ_yield_stress:", integ_yield_stress)

    element_ctr = zeros(nElement,3);
    #element_ctr = cal_element_ctr(coordmat, elementmat);

    output_num = 100
    d_out = Int(floor(time_num / output_num))

    output_disp = Array{Float64,2}(undef, fn, output_num)
    output_element_flag = Array{Int,2}(undef, nElement, output_num)
    output_data = Array{IntegDataType,1}(undef,output_num)

    integ_data = IntegDataType(0, integ_stress, integ_strain, integ_plastic_strain, integ_eq_plastic_strain, integ_triax_stress)
    node_data = cal_node_stress_strain(nNode, elementmat, integ_num, integ_data)
    write_vtk(0, coordmat, elementmat, element_flag, disp, velo, node_data)


    output = [];

    println("")
    i_out = 1
@time for t = 1 : time_num

        if rem(t,100) == 0
          #println(t)
          str_ = @sprintf("\r%.4e / %.4e     ", t*d_time, end_time)
          print(str_)
        end

        #set force
        #external_force = zeros(fn, 1);
        external_force .= 0.0      # -> 0 alloc
        external_force[ force_dof ] = force_v;

        if contact_flag >= 1
            c_force3 .= 0.0  #->0 alloc
            #d_node .= 0.0 #->0 alloc

            if gpu_contact_flag == 0
                #@time 
                cal_contact_force(c_force3, CT, instance_pair, cp_index, #CV, 
                                position, velo, diag_M, elementMinSize, elementMaxSize, d_max,  # d_node, d_node_pre, 
                                element_flag, elementmat, bug_report, t*d_time);
                # -> 0 alloc

                if Nth > 1
                    @inbounds @floop for  i = 1 : fn
                        @inbounds for  j = 2 : Nth
                            c_force3[i,1] += c_force3[i,j] 
                        end
                    end
                end

            else

                cal_contact_force_gpu(c_force3, CT, instance_pair, cp_index, #CV, 
                                position, velo, diag_M, elementMinSize, elementMaxSize, d_max,  
                                element_flag, elementmat, bug_report, t*d_time);
                #println("c_force3:", c_force3)
            end

        #=    print(bug_report, @sprintf("time=%.16e\n", t*d_time) ) 
            @inbounds for  i = 1 : fn
                if c_force3[i] != 0.0
                    print(bug_report, @sprintf("%.16e\n", c_force3[i]) ) 
                end
            end
        =#

                #external_force .+= c_force3
            @inbounds @floop for  i = 1 : fn
                external_force[i] += c_force3[i]
            end

            fmax, idx_fmax = findmax(abs.(external_force))
            #print(bug_report, @sprintf("time=%.6e, fmax=%.6e, index=%d\n", 
            #                            t*d_time, fmax, idx_fmax) )


            #@inbounds @floop for  i = 1 : nNode  
            #    d_node_pre[i] = d_node[i]
            #end
        end
  

        #update position
            #disp_new = 1.0 ./ (diag_M/d_time^2 .+ diag_C/2.0/d_time) .* ( external_force .- Q .+ diag_M/d_time^2.0 .* (2.0*disp .- disp_pre) .+ diag_C/2.0/d_time .* disp_pre );
        #disp_new .= 1.0 ./ (diag_M/d_time^2 .+ diag_C/2.0/d_time) .* ( external_force .- Q .+ diag_M/d_time^2.0 .* (2.0*disp .- disp_pre) .+ diag_C/2.0/d_time .* disp_pre );
        # .=  -> much faster than = 
        # But, 140k alloc

            #println("disp_new:")
            #println(disp_new)

        #println("disp_new_i:")
        # @floop 
        @inbounds @floop for  i = 1 : fn    # much faster than [ =  dot operations ],     Almost equal to .=
                #disp_new[i] = 1.0 / (diag_M[i]*div_d_time2 + diag_C[i]/2.0*div_d_time) * ( external_force[i] - Q[i] + diag_M[i]*div_d_time2 * (2.0*disp[i] - disp_pre[i]) + diag_C[i]/2.0*div_d_time * disp_pre[i] );
            disp_new[i] = 1.0 / (diag_M[i]/d_time^2 + diag_C[i]/2.0/d_time) * ( external_force[i] - Q[i] + diag_M[i]/d_time^2.0 * (2.0*disp[i] - disp_pre[i]) + diag_C[i]/2.0/d_time * disp_pre[i] );
                #print(v,", ")
                #disp_new[i] = v
        end  # -> 0 alloc
        #println("")
        #println(disp_new)

#=        println("fn:", fn)
        println("disp_new:", size(disp_new))
        println("diag_M:", size(diag_M))
        println("diag_C:", size(diag_C))
        println("external_force:", size(external_force))
        println("Q:", size(Q))
        println("disp:", size(disp))
        println("disp_pre:", size(disp_pre)) =#

        #println("")

        #% Boundary Conditions
        #%disp_new( fix_dof ) = 0;
        #%disp_new( input_disp_dof ) = input_disp_mag * t / time_num ;
        for i = 1 : length(MODEL.BC)   # -> 0 alloc
            amp = 1.0;
            if length(MODEL.BC[i].amp_name) > 0  #length(MODEL.BC[i].amplitude) > 0
                time_index = 1;
                current_time = t * d_time;
                a_t = MODEL.BC[i].amplitude.time;
                a_v = MODEL.BC[i].amplitude.value;
                for j = 1 : length(a_t) - 1
                    if current_time >= a_t[j] && current_time <= a_t[j+1]
                        time_index = j;
                        break
                    end
                end
                amp = a_v[time_index] + (a_v[time_index+1] - a_v[time_index]) * 
                         (current_time - a_t[time_index]) / (a_t[time_index+1]-a_t[time_index]);
            end

            #println("BC[", i, "]:", MODEL.BC[i])
            #println("size(MODEL.BC[i].dof):", size(MODEL.BC[i].dof) )

            for j = 1 : length(MODEL.BC[i].dof)
                dof = MODEL.BC[i].dof[j]  #::Vector{Int};
                v = MODEL.BC[i].value[j];
                #println("dof:", dof)
                #println("v:", v)
                #println("amp:", amp)
                
                disp_new[ dof ] .= v * amp;
                #for k = 1 : length(dof)   # -> slow
                #    disp_new[ dof[k] ] = v * amp;
                #end
            end
        end

              #d_disp .= copy(disp_new) .- copy(disp);  # -> 40k alloc
              #disp_pre = copy(disp);  # -> 20k alloc
              #disp = copy(disp_new);  # -> 20k alloc
              #velo .= copy(d_disp) ./ d_time;  # -> 20k alloc

        @inbounds @floop for  i = 1 : fn  # -> 0 alloc
            d_disp[i] = disp_new[i] - disp[i]
            disp_pre[i] = disp[i]
            disp[i] = disp_new[i]
            velo[i] = d_disp[i] / d_time
        end


#=          mean_velo = sum(abs.(velo)) / fn;
            for i = 1 : fn
                if abs(velo[i]) < mean_velo * 1.0E-03
                 velo[i] = 0.0;
                end
            end =#

            #d_disp3 = reshape(d_disp,3,nNode)  #-> 40k alloc
            #for i = 1 : nNode # -> 480k alloc
            #    d_disp_norm[i] = my3norm(d_disp3[:,i])
            #end

        @inbounds @floop for i = 1 : nNode  #-> 0 alloc
            x = d_disp[i*3-2]
            y = d_disp[i*3-1]
            z = d_disp[i*3]
            d_disp_norm[i] = sqrt(x*x+y*y+z*z)

            position[1,i] = coordmat[1,i] + disp[i*3-2]
            position[2,i] = coordmat[2,i] + disp[i*3-1]
            position[3,i] = coordmat[3,i] + disp[i*3]
        end
        
        #d_max = maximum(d_disp_norm)
            #println("d_max:", d_max)
        d_max, idx_dmax = findmax(d_disp_norm)
            #print(bug_report, @sprintf("time=%.6e, d_max=%.6e, index=%d, coordmat:%.6e,%.6e,%.6e\n", 
            #                            t*d_time, d_max, idx_dmax, coordmat[1,idx_dmax],coordmat[2,idx_dmax],coordmat[3,idx_dmax]) ) 
        

        Qe .= 0.0
        #@time 
        cal_stress_hexa(#Barray, BVarray, #detJarray, #e_position,  #Dmat, 
                        Qe, integ_stress, integ_strain, integ_yield_stress, integ_eq_plastic_strain, 
                        position, d_disp, elementmat, element_flag, integ_num, Pusai_mat,  
                        MODEL.MATERIAL, element_material, elementMinSize, elementVolume); # ->OK
        Q .= 0.0                        
        @inbounds for e = 1 : nElement
            @inbounds for i = 1 : 8
                Q[ 1 + (elementmat[i,e]-1)*3 ] += Qe[1 + (i-1)*3, e];
                Q[ 2 + (elementmat[i,e]-1)*3 ] += Qe[2 + (i-1)*3, e];
                Q[ 3 + (elementmat[i,e]-1)*3 ] += Qe[3 + (i-1)*3, e];
            end
        end
  
        cal_triax_stress( integ_stress, integ_triax_stress); #-> less than 10% of whole
        #-> 0 alloc

        #print(bug_report, @sprintf("time=%.6e, Vemin=%.6e\n", t*d_time, minimum(elementVolume)) ) 

        deleted_element = Int[];
        
        if flag_fracture == 1  #->20k

        #=          %--- fracture stress ---%
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
                    % end =#
 
            #%--- fracture strain ---%
            for i = 1 : nElement
              mat_id = element_material[i];
              
              nd =size(MODEL.MATERIAL[mat_id].ductile, 1);
              if nd > 0
  
                  #v_e = sum( integ_eq_plastic_strain[(1:integ_num).+(i-1)*integ_num] ) / integ_num;
                  #t_e = sum( integ_triax_stress[(1:integ_num).+(i-1)*integ_num] ) / integ_num;
                  # -> huge alloc

                  v_e = 0.0
                  t_e = 0.0
                  @inbounds for j = 1 : integ_num # 8
                        v_e += integ_eq_plastic_strain[j+(i-1)*integ_num]
                        t_e += integ_triax_stress[j+(i-1)*integ_num]
                  end
                  v_e /= integ_num
                  t_e /= integ_num

                  if t_e < 0
                      continue
                  end
    
                  ductile_ = MODEL.MATERIAL[mat_id].ductile;
                  fr_e = ductile_[nd,1];
                  for j = 1 : nd-1
                      if t_e >= ductile_[j,2] && t_e < ductile_[j+1, 2]
                          fr_e = ductile_[j,1] + (ductile_[j+1,1] - ductile_[j,1]) / (ductile_[j+1,2] - ductile_[j,2]) * (t_e - ductile_[j,2]);
                          break
                      end
                  end
 
                  if v_e >= fr_e  &&  element_flag[i] == 1
                      element_flag[i] = 0;
                      append!(deleted_element, i)   # -> 1 alloc
                      println("Element deleted:", sum(element_flag), "/", nElement)
                            #[sum(element_flag)  nElement]
                            #clstr_num = 0;

                      #integ_stress[:, (1:integ_num).+(i-1)*integ_num] .= 0;
                      #integ_strain[:, (1:integ_num).+(i-1)*integ_num] .= 0;
                      @inbounds for j = 1 : integ_num
                            integ_stress[1, j+(i-1)*integ_num] = 0.0
                            integ_stress[2, j+(i-1)*integ_num] = 0.0
                            integ_stress[3, j+(i-1)*integ_num] = 0.0
                            integ_stress[4, j+(i-1)*integ_num] = 0.0
                            integ_stress[5, j+(i-1)*integ_num] = 0.0
                            integ_stress[6, j+(i-1)*integ_num] = 0.0

                            integ_strain[1, j+(i-1)*integ_num] = 0.0
                            integ_strain[2, j+(i-1)*integ_num] = 0.0
                            integ_strain[3, j+(i-1)*integ_num] = 0.0
                            integ_strain[4, j+(i-1)*integ_num] = 0.0
                            integ_strain[5, j+(i-1)*integ_num] = 0.0
                            integ_strain[6, j+(i-1)*integ_num] = 0.0
                      end
                        
                  end
  
              end
  
            end
        
        end
  
        #%--- Update surface ---%
        if contact_flag > 0
            for i = deleted_element
    
                instance_id = element_instance[i]
                ele_id = i - MODEL.INSTANCE[instance_id].element_offset
                add_c_triangles, add_c_triangles_eleid, add_c_nodes = add_surface_triangle(MODEL.INSTANCE[instance_id], ele_id)
                #println("add_c_triangles:", add_c_triangles)
                #println("add_c_triangles_eleid:", add_c_triangles_eleid)
                #println("add_c_nodes:", add_c_nodes)
    
                for c = 1 : length(instance_pair)
    
                    #% j -> triangle,   i -> point
                    i_instance = instance_pair[c][1]
                    j_instance = instance_pair[c][2]
                    #println("i_instance:", i_instance, " j_instance:", j_instance)
                        
                    if i_instance == instance_id  #i_instance
                        #CT[c].c_nodes_i = c_nodes .+ MODEL.INSTANCE[i].node_offset  
                        append!(CT[c].c_nodes_i, add_c_nodes .+ MODEL.INSTANCE[instance_id].node_offset )
                        unique!(CT[c].c_nodes_i)
                        
                    elseif j_instance == instance_id  
                        #CT[c].c_nodes_j = c_nodes .+ MODEL.INSTANCE[i].node_offset
                        #CT[c].c_triangles = c_triangles .+ MODEL.INSTANCE[i].node_offset
                        #CT[c].c_triangles_eleid = c_triangles_eleid .+ MODEL.INSTANCE[i].element_offset
    
                        append!(CT[c].c_nodes_j, add_c_nodes .+ MODEL.INSTANCE[instance_id].node_offset )
                        unique!(CT[c].c_nodes_j)
                        append!(CT[c].c_triangles_eleid, add_c_triangles_eleid .+ MODEL.INSTANCE[instance_id].element_offset )
                        CT[c].c_triangles = vcat(CT[c].c_triangles, add_c_triangles .+ MODEL.INSTANCE[instance_id].node_offset)
                    end
         
                    #println("CT[", c, "]:", CT[c])  # -> Match with Matlab
                end
    
            end
        end

#=        
        updated_instance = element_instance[deleted_element]  # -> 1 alloc
        unique!(updated_instance) # -> 0 alloc
        
        if length(updated_instance) > 0  && contact_flag == 1   # -> small time cost
            #println("updated_instance")
            for i = updated_instance
  
#=                  #%--- INSTANCE base 
                    u_ele = Int[];
                    for k = 1 : MODEL.INSTANCE[i].nElement
                        offset = MODEL.INSTANCE[i].element_offset;
                        if element_flag[k+offset] == 1
                            append!(u_ele, k)
                        end
                    end
        
                    #c_triangles, c_triangles_eleid, c_nodes = get_surface_triangle(MODEL, i, u_ele);
                        #MODEL.INSTANCE[i].c_triangles = c_triangles;
                        #MODEL.INSTANCE[i].c_triangles_eleid = c_triangles_eleid;
                        #MODEL.INSTANCE[i].c_nodes = c_nodes;
                    =#
                
                
                u_ele = Int[];
                for k = 1 : MODEL.INSTANCE[i].nElement 
                   offset = MODEL.INSTANCE[i].element_offset; 
                   if element_flag[k+offset] == 1
                        append!(u_ele, k)
                   end
                end

                if all_exterior_flag == 1
                    
                    c_triangles, c_triangles_eleid, c_nodes = get_surface_triangle(MODEL.INSTANCE[i], u_ele, MODEL.INSTANCE[i].elements)
                    
                    for c = 1 : length(instance_pair)

                        #cc = cp_index[c]
                
                        #% j -> triangle,   i -> point
                        i_instance = instance_pair[c][1]
                        j_instance = instance_pair[c][2]
                        #println("i_instance:", i_instance, " j_instance:", j_instance)
                            
                        if i_instance == i  #i_instance
                            CT[c].c_nodes_i = c_nodes .+ MODEL.INSTANCE[i].node_offset                            
                        elseif j_instance == i  
                            CT[c].c_nodes_j = c_nodes .+ MODEL.INSTANCE[i].node_offset
                            CT[c].c_triangles = c_triangles .+ MODEL.INSTANCE[i].node_offset
                            CT[c].c_triangles_eleid = c_triangles_eleid .+ MODEL.INSTANCE[i].element_offset
                        end
            
                        #CT[c].c_nodes_i = c_nodes_i
                        #CT[c].c_nodes_j = c_nodes_j
                        #CT[c].c_triangles = c_triangles
                        #CT[c].c_triangles_eleid = c_triangles_eleid
                        
                        #println("CT[", c, "]:", CT[c])  # -> Match with Matlab
                    end

                end # all_exterior_flag == 1

                if all_exterior_flag == 0
                    # if define local contact, no update for contact surface triangle

#=                  for c = 1 : length(cp_index)

                        cc = cp_index[c]
                
                        #% j -> triangle,   i -> point
                        #i_instance = instance_pair[c][1]
                        #j_instance = instance_pair[c][2]
                        #println("i_instance:", i_instance, " j_instance:", j_instance)
                            
                        if MODEL.CP[cc].instance_id_1 == i  #i_instance
                            c_triangles, c_triangles_eleid, c_nodes = get_surface_triangle(MODEL.INSTANCE[i], u_ele, MODEL.CP[cc].elements_1);
                            CT[c].c_nodes_i = c_nodes .+ MODEL.INSTANCE[i].node_offset                            
                        elseif MODEL.CP[cc].instance_id_2 == i  
                            c_triangles, c_triangles_eleid, c_nodes = get_surface_triangle(MODEL.INSTANCE[i], u_ele, MODEL.CP[cc].elements_2);
                            CT[c].c_nodes_j = c_nodes .+ MODEL.INSTANCE[i].node_offset
                            CT[c].c_triangles = c_triangles .+ MODEL.INSTANCE[i].node_offset
                            CT[c].c_triangles_eleid = c_triangles_eleid .+ MODEL.INSTANCE[i].element_offset
                        end
            
                        #CT[c].c_nodes_i = c_nodes_i
                        #CT[c].c_nodes_j = c_nodes_j
                        #CT[c].c_triangles = c_triangles
                        #CT[c].c_triangles_eleid = c_triangles_eleid
                        
                        #println("CT[", c, "]:", CT[c])  # -> Match with Matlab
                    end
=#
                end # all_exterior_flag == 0

                #%--- Contact Pair base                
#=                  for j = 1 : length( MODEL.CP )
    
                        if MODEL.CP[j].instance_id_1 == i          
                            c_triangles, c_triangles_eleid, c_nodes = get_surface_triangle(MODEL.INSTANCE[i], u_ele, MODEL.CP[j].elements_1);
                            MODEL.CP[j].c_triangles_1 = c_triangles;
                            MODEL.CP[j].c_triangles_eleid_1 = c_triangles_eleid;
                            MODEL.CP[j].c_nodes_1 = c_nodes;
                            #println("c_triangles_1:", size(c_triangles,1),":",c_triangles)
                            #println("c_triangles_eleid_1:", length(c_triangles_eleid),":",c_triangles_eleid)
                            #println("c_nodes_1:", length(c_nodes),":",c_nodes)
                            
                        end
    
                        if MODEL.CP[j].instance_id_2 == i            
                            c_triangles, c_triangles_eleid, c_nodes = get_surface_triangle(MODEL.INSTANCE[i], u_ele, MODEL.CP[j].elements_2);
                            MODEL.CP[j].c_triangles_2 = c_triangles;
                            MODEL.CP[j].c_triangles_eleid_2 = c_triangles_eleid;
                            MODEL.CP[j].c_nodes_2 = c_nodes;
                            #println("c_triangles_2:", size(c_triangles,1),":",c_triangles)
                            #println("c_triangles_eleid_2:", length(c_triangles_eleid),":",c_triangles_eleid)
                            #println("c_nodes_2:", length(c_nodes),":",c_nodes)
                            # -> match with Matlab
                        end

                    end #for j
=#
            end # for i
        end # if
=#  

        if rem(t,d_out) == 0   
            output_disp[:, i_out] = disp  # -> 0 alloc
            output_element_flag[:, i_out] = element_flag  # -> 0 alloc
            output_data[i_out] = IntegDataType(i_out, integ_stress, integ_strain, integ_plastic_strain, integ_eq_plastic_strain, integ_triax_stress)
            # -> 1 alloc
            
            node_value = cal_node_stress_strain(nNode, elementmat, integ_num, output_data[i_out] );  #14.4k / 40 alloc
            write_vtk(i_out, coordmat, elementmat, output_element_flag[:,i_out], output_disp[:,i_out], copy(velo), node_value); #  41.28k / 40 alloc

            i_out = i_out + 1;
        end

        #position
        #output = push!(output, disp[fn]);


        if t == time_num
            println("")
        end
    end
    println("")
    
    #println("position:", position)
    #println("disp3:", reshape(disp,3,nNode)') # ->OK
    #println("Q:", Q)
    #println("integ_triax_stress:", integ_triax_stress)
    #println("c_force3:", c_force3)
    #println("external_force:", external_force)
    #println("d_node:", d_node)

    #drawElement(coordmat,elementmat)
    #drawGraph(output)

    #println("integ_stress:", integ_stress) # ->OK
    #println("node_stress:", node_stress) # ->OK
    #println("node_mises_stress:", node_mises_stress)

#=    for i = 1 : output_num
        node_value = cal_node_stress_strain(nNode, elementmat, integ_num, output_data[i] );
        write_vtk(i, coordmat, elementmat, output_element_flag[:,i], output_disp[:,i], node_value);
    end=#

    close(bug_report)

    return output

end


#function cal_triax_stress( integ_stress::Array{Float64,2} )
function cal_triax_stress( integ_stress::Array{Float64,2}, integ_triax_stress::Array{Float64,1} )
    #n = size(integ_stress,1);
    n = size(integ_stress,2);
    #integ_triax_stress = zeros(n);
    integ_triax_stress .= 0.0

    #T = @SMatrix zeros(3,3)
    #p = @SVector zeros(3) #Array{Float64}(undef,3)
    #ox = 0.0
    #p = @MVector zeros(3) # -> no meaning

    #div3::Float64 = 1.0/3.0

    @inbounds @floop for i =  1 : n
        ox = integ_stress[1,i]
        oy = integ_stress[2,i]
        oz = integ_stress[3,i]
        txy = integ_stress[4,i]
        tyz = integ_stress[5,i]
        txz = integ_stress[6,i]
        #%oeq = sqrt(0.5*( (ox-oy)^2 + (oy-oz)^2 + (ox-oz)^2 + 6*(txy^2 + tyz^2 + txz^2) ));

        T = @SMatrix  [ox  txy txz
                       txy oy  tyz
                       txz tyz oz ];
        p = eigvals(T) 
        #%v = (p(1)+p(2)+p(3))/3 / oeq;

        oeq = sqrt(0.5*( (p[1]-p[2])^2 + (p[2]-p[3])^2 + (p[3]-p[1])^2  ));
        #oeq = sqrt(0.5*( (p[1]-p[2])*(p[1]-p[2]) + (p[2]-p[3])*(p[2]-p[3]) + (p[3]-p[1])*(p[3]-p[1])  ));
        if oeq < 1E-10  #%oeq == 0
            continue
        end

        v = (p[1]+p[2]+p[3]) / 3.0 / oeq;   # /3
        integ_triax_stress[i] = v;
    end

    #return integ_triax_stress::Array{Float64,1}

end


#function cal_stress_hexa(position_, d_disp_, elementmat, element_flag, integ_num, Pusai_mat, integ_stress_pre_, 
#                         MATERIAL, element_material, integ_yield_stress_, integ_eq_plastic_strain_, elementMinSize)
#function cal_stress_hexa(Q, integ_strain_,
#                         position_, d_disp_, elementmat, element_flag, integ_num, Pusai_mat, integ_stress_pre_, 
#                         MATERIAL, element_material, integ_yield_stress_, integ_eq_plastic_strain_, elementMinSize)
#function cal_stress_hexa(Q, integ_stress_pre_, integ_strain_, integ_yield_stress_, integ_eq_plastic_strain_, 
#                         position_, d_disp_, elementmat, element_flag, integ_num, Pusai_mat,  
#                         MATERIAL, element_material, elementMinSize)
function cal_stress_hexa(#Barray, BVarray, #detJarray, #e_position,  #Dmat,   #-> huge reduce alloc and memory
                         Qe, integ_stress_pre_, integ_strain_, integ_yield_stress_, integ_eq_plastic_strain_, 
                         position_, d_disp_, elementmat, element_flag, integ_num, Pusai_mat,  
                         MATERIAL, element_material, elementMinSize, elementVolume)

    nElement = length(element_flag);  #size(elementmat,1);
            #d_integ_stress_ = zeros(6, nElement*integ_num);
            #d_integ_strain_ = zeros(6, nElement*integ_num);
            #d_integ_yield_stress_ = zeros(nElement*integ_num);
            #d_integ_plastic_strain_ = zeros(6, nElement*integ_num);
            #d_integ_eq_plastic_strain_ = zeros(nElement*integ_num);

    fn = size(position_,1)*3;
            #Q = zeros(fn)
    #Q .= 0.0;

            # weight = 1
            # W = [2.0]; # integ_num = 1
            #if integ_num == 8
            #    W = @SVector ones(8);
            #end
    W = 1.0
    #div3::Float64 = 1.0/3.0

        #Barray = zeros(6,24,integ_num);  #@MArray -> huge increase alloc and memory
        #BVarray = zeros(6,24,integ_num);
    
    #Barray = [ @MMatrix zeros(6,24) for i=1:8] # -> huge reduce alloc and memory compared to zeros(6,24,integ_num)
    #BVarray = [ @MMatrix zeros(6,24) for i=1:8]
    #for i = 1 : 8
    #    Barray[i] .= 0.0
    #    BVarray[i] .= 0.0
    #end

        #d_u = @MVector zeros(24); # Memory reduce
        #B = @MMatrix zeros(6,24);
        #P2 = @MMatrix zeros(3,8);
            #Jarray = zeros(3,3,integ_num);
    
        #BVbar = @MMatrix zeros(6,24); # @MMatrix -> reduce alloc and memory
        #BVbar_i = @MMatrix zeros(6,24);
        #BV_i = @MMatrix zeros(6,24);
        #Bfinal = @MMatrix zeros(6,24); # @MMatrix -> huge reduce alloc and memory
        #q_vec = @MVector zeros(24); # Memory reduce
        #q_vec_i = @MVector zeros(24)
        #tri_dev_stress = @MVector zeros(6); #@MVector -> huge reduce alloc and memory
        #pre_stress = @MVector zeros(6); #@MVector -> bit reduce alloc and memory
        #tri_stress = @MVector zeros(6); #@MVector -> reduce alloc and memory
        #mean_stress_vec = @MVector zeros(6); #@MVector -> huge reduce alloc and memory
        #final_dev_stress = @MVector zeros(6);
        #final_stress = @MVector zeros(6); #@MVector -> bit reduce alloc and memory
        #d_o_vec = @MVector zeros(6); #@MVector -> bit reduce alloc and memory
        #d_e_vec = @MVector zeros(6); #@MVector -> even
            #o_vec = zeros(6); #@MVector -> huge increase alloc and memory
            #e_position = zeros(8,3); # @MMatrix -> huge increase alloc and memory, causes compilation time


            #d_e_vec = @MVector zeros(Float64, 6) # -> no meaning
            #J = zeros(3,3)

            # cause small alloc
                #detJarray = zeros(integ_num);  #@MVector -> huge increase alloc and memory
                #e_position = zeros(3,8)  # -> 1 alloc
                #Dmat = @MMatrix zeros(6,6)  #@MMatrix -> huge reduce alloc and memory,  # -> 1 alloc
                #Dmat .= MATERIAL[mat_id].Dmat;

    #=
    mat_id = element_material[1]
    G = MATERIAL[mat_id].G;
    Hd = MATERIAL[mat_id].Hd
    plastic_property_ = MATERIAL[mat_id].plastic;
    npp = size(plastic_property_,1)

    @inbounds for i = 1 : 36  # -> 0 alloc
        Dmat[i] = MATERIAL[mat_id].Dmat[i]
    end
    =#

    #q_mat = @MMatrix zeros(24, nElement)
    
    #  
    @inbounds @floop for e = 1 : nElement  #

        if element_flag[e] == 0
            continue;
        end

        mat_id = element_material[e];
            #Dmat = MATERIAL[mat_id].Dmat;
        G = MATERIAL[mat_id].G;
        plastic_property_ = MATERIAL[mat_id].plastic;
        Hd = MATERIAL[mat_id].Hd
        npp = size(plastic_property_,1)
        Dmat =  @MMatrix zeros(6,6)
        @inbounds for i = 1 : 36  # -> 0 alloc
            Dmat[i] = MATERIAL[mat_id].Dmat[i]
        end

    #=      if mat_id != element_material[e]
                mat_id = element_material[e]
                    #Dmat .= MATERIAL[mat_id].Dmat;
                @inbounds for i = 1 : 36  # -> 0 alloc
                    Dmat[i] = MATERIAL[mat_id].Dmat[i]
                end
                G = MATERIAL[mat_id].G;
                plastic_property_ = MATERIAL[mat_id].plastic;
                Hd = MATERIAL[mat_id].Hd
                npp = size(plastic_property_,1)
            end=#

                #e_position = position_[ elementmat[e,:], :];
                #println("e_position:", e_position)

            #d_u = zeros(24,1);
        d_u = @MVector zeros(24)
        e_position = @MMatrix  zeros(3,8)  #
        @inbounds for i = 1 : 8   # 0 alloc
            d_u[1+(i-1)*3] = d_disp_[ 1+(elementmat[i,e]-1)*3 ];  # -> 0 alloc
            d_u[2+(i-1)*3] = d_disp_[ 2+(elementmat[i,e]-1)*3 ];
            d_u[3+(i-1)*3] = d_disp_[ 3+(elementmat[i,e]-1)*3 ];

            e_position[1,i] = position_[ 1, elementmat[i,e]];  # -> reduce of alloc and memory, comparing "e_position = position_[ elementmat[e,:], :]"
            e_position[2,i] = position_[ 2, elementmat[i,e]];
            e_position[3,i] = position_[ 3, elementmat[i,e]];
        end

                #if norm(d_u) < elementMinSize * 1.0E-7
                #    continue
                #end

        q_vec = @MVector zeros(24)
            #BVbar .= 0.0
        BVbar = @MMatrix zeros(6,24)
            #V = 0.0;

        V = cal_BVbar_hexa( Pusai_mat, e_position, BVbar)   
        elementVolume[e] = V

#=        Parray = [ @MMatrix zeros(3,8) for i=1:8]
        detJarray = @MVector zeros(8)
        @inbounds for i = 1 : integ_num

                    #detJi = cal_B_BVbar_hexa( Pusai_mat[i], e_position, B, BV_i, BVbar_i)   
                    #Barray[i] .= B
                    #BVarray[i] .= BV_i
                    
                    #detJi = cal_B_BVbar_hexa( Pusai_mat[i], e_position, Barray[i], BVarray[i], BVbar_i)   
                    #BVbar .= BVbar .+ BVbar_i
                #detJi = cal_B_BVbar_hexa( Pusai_mat[i], e_position, Barray[i], BVarray[i], BVbar)   
                #V = V + detJi
            #detJi = cal_B_BVbar_hexa( Pusai_mat[i], e_position, Barray[i], BVarray[i], Parray[i])   
            detJi = cal_P_hexa( Pusai_mat[i], e_position, Parray[i])   
            detJarray[i] = detJi
            
        end =#
        #BVbar .= BVbar / V;

        
        Bfinal = @MMatrix zeros(6,24)

        @inbounds for i = 1 : integ_num

                #Bfinal .= Barray[i] .+ BVbar .- BVarray[i];
            #Bfinal = Barray[i] .+ BVbar .- BVarray[i];
          #Bfinal = @MMatrix zeros(6,24)
          Bfinal .= 0
          #cal_Bfinal(Bfinal, Parray[i], BVbar)
          detJ = cal_Bfinal(Bfinal, BVbar, Pusai_mat[i], e_position)

            #d_e_vec .= Bfinal * d_u;  # .=  -> huge reduce for alloc and memory
            #d_o_vec .= Dmat * d_e_vec;
          d_e_vec = Bfinal * d_u
          d_o_vec = Dmat * d_e_vec;

          index_i = (e-1)*integ_num+i;
            #pre_stress = integ_stress_pre_[:,index_i];
            #pre_stress = copy(integ_stress_[:,index_i]);
            #pre_stress = integ_stress_pre_[index_i,:];
          pre_stress = @MVector zeros(6)
          pre_stress[1] = integ_stress_pre_[1,index_i];  # -> reduce of alloc and memory, comparing "pre_stress = integ_stress_pre_[:,index_i]"
          pre_stress[2] = integ_stress_pre_[2,index_i];
          pre_stress[3] = integ_stress_pre_[3,index_i];
          pre_stress[4] = integ_stress_pre_[4,index_i];
          pre_stress[5] = integ_stress_pre_[5,index_i];
          pre_stress[6] = integ_stress_pre_[6,index_i]; 

            #final_stress .= pre_stress .+ d_o_vec
          final_stress = pre_stress .+ d_o_vec

          #println("d_e_vec:", d_e_vec)
          #println("d_o_vec:", d_o_vec)
          #println("pre_stress:", pre_stress)
          #println("pre_stress':", pre_stress')

          if length( plastic_property_ ) > 0 
                #tri_stress .= pre_stress .+ d_o_vec;  # .= .+  -> huge reduce of alloc and memory
               tri_stress = pre_stress .+ d_o_vec;
               mean_stress = ( tri_stress[1]+tri_stress[2]+tri_stress[3] ) / 3.0;
                # @SVector # -> even
                #=tri_dev_stress = [tri_stress[1] - mean_stress
                                    tri_stress[2] - mean_stress
                                    tri_stress[3] - mean_stress
                                    tri_stress[4]
                                    tri_stress[5]
                                    tri_stress[6]];=#
               tri_dev_stress = @MVector zeros(6)
               tri_dev_stress[1] = tri_stress[1] - mean_stress; # even
               tri_dev_stress[2] = tri_stress[2] - mean_stress;
               tri_dev_stress[3] = tri_stress[3] - mean_stress;
               tri_dev_stress[4] = tri_stress[4];
               tri_dev_stress[5] = tri_stress[5];
               tri_dev_stress[6] = tri_stress[6];

               mean_stress_vec = @MVector zeros(6)
               mean_stress_vec[1] = mean_stress
               mean_stress_vec[2] = mean_stress
               mean_stress_vec[3] = mean_stress

               tri_mises_stress = sqrt( 1.5 * (tri_dev_stress[1]^2 + tri_dev_stress[2]^2 + tri_dev_stress[3]^2 + 
                                               2*tri_dev_stress[4]^2 + 2*tri_dev_stress[5]^2 + 2*tri_dev_stress[6]^2) );
               y = integ_yield_stress_[ index_i ];
               if tri_mises_stress > y
                   p_index = 1;
                   for j = 2 : npp # size(plastic_property_,1)
                       if integ_eq_plastic_strain_[ index_i ] <= plastic_property_[j,2]
                            p_index = j-1;
                            break
                       end
                       if j == npp # size(plastic_property_,1)
                            p_index = j-1;
                       end
                   end
                    #%p_index
                    #H = (plastic_property_[p_index+1,1] - plastic_property_[p_index,1]) / (plastic_property_[p_index+1,2] - plastic_property_[p_index,2]);
                   H = Hd[p_index]
                   d_ep = ( tri_mises_stress - y ) / (3*G + H);
                   
                    #d_integ_eq_plastic_strain_[ index_i ] = d_ep;
                    #d_integ_yield_stress_[ index_i ] = H * d_ep;
                   
                    #final_dev_stress .= tri_dev_stress * (y+H*d_ep) / tri_mises_stress;
                   final_dev_stress = tri_dev_stress * (y+H*d_ep) / tri_mises_stress;
                    #final_stress = final_dev_stress + @SVector [mean_stress, mean_stress, mean_stress, 0, 0, 0];
                    #final_stress .= final_dev_stress .+ [mean_stress, mean_stress, mean_stress, 0, 0, 0];
                   final_stress .= final_dev_stress .+ mean_stress_vec;
                   
                   #d_o_vec .= final_stress .- pre_stress; # .= .-  -> even

                   integ_eq_plastic_strain_[ index_i ] += d_ep;
                   integ_yield_stress_[ index_i ] += H * d_ep;
               #else
                   #final_stress .= pre_stress .+ d_o_vec
               end # if

          #else
                #final_stress .= pre_stress .+ d_o_vec
          end #if

                #d_integ_stress_[index_i,:] = d_o_vec'
                #d_integ_strain_[index_i,:] = d_e_vec'
                #d_integ_stress_[:,index_i] = d_o_vec
                #d_integ_strain_[:,index_i] = d_e_vec
                #d_integ_strain_[1,index_i] = d_e_vec[1]  # -> same
                #d_integ_strain_[2,index_i] = d_e_vec[2]
                #d_integ_strain_[3,index_i] = d_e_vec[3]
                #d_integ_strain_[4,index_i] = d_e_vec[4]
                #d_integ_strain_[5,index_i] = d_e_vec[5]
                #d_integ_strain_[6,index_i] = d_e_vec[6]

                #integ_stress_[:,index_i] += d_o_vec
                #integ_strain_[:,index_i] += d_e_vec  # -> huge allocations
          integ_strain_[1,index_i] += d_e_vec[1]
          integ_strain_[2,index_i] += d_e_vec[2]
          integ_strain_[3,index_i] += d_e_vec[3]
          integ_strain_[4,index_i] += d_e_vec[4]
          integ_strain_[5,index_i] += d_e_vec[5]
          integ_strain_[6,index_i] += d_e_vec[6]

        #=        integ_stress_pre_[1,index_i] += d_o_vec[1]
                integ_stress_pre_[2,index_i] += d_o_vec[2]
                integ_stress_pre_[3,index_i] += d_o_vec[3]
                integ_stress_pre_[4,index_i] += d_o_vec[4]
                integ_stress_pre_[5,index_i] += d_o_vec[5]
                integ_stress_pre_[6,index_i] += d_o_vec[6] =#

          integ_stress_pre_[1,index_i] = final_stress[1]
          integ_stress_pre_[2,index_i] = final_stress[2]
          integ_stress_pre_[3,index_i] = final_stress[3]
          integ_stress_pre_[4,index_i] = final_stress[4]
          integ_stress_pre_[5,index_i] = final_stress[5]
          integ_stress_pre_[6,index_i] = final_stress[6] 

                #o_vec .= pre_stress .+ d_o_vec;  # .= .+  -> reduce of alloc and memory

                #q_vec_i = B' * o_vec;
                #q_vec_i = Bfinal' * o_vec; # .= -> huge increase alloc and memory
            #q_vec_i .= Bfinal' * final_stress # .=   -> huge reduce of alloc and memory
          q_vec_i = Bfinal' * final_stress

                #q_vec = q_vec + W[i]*W[i]*W[i] * detJarray[i] * q_vec_i;
                #q_vec .+=  W[i]*W[i]*W[i] * detJarray[i] * q_vec_i; # .= .+  -> huge reduce of alloc and memory  #->0.8M alloc
          @inbounds for j = 1 : 24   # -> 1 alloc -> 0 alloc
                Qe[j, e] += W * W * W * detJ * q_vec_i[j]
                #q_vec[j] += W * W * W * detJ * q_vec_i[j]
                    #q_vec[j] += W * W * W * detJarray[i] * q_vec_i[j]
                            #q_vec[j] += W[i]*W[i]*W[i] * detJarray[i] * q_vec_i[j] # -> 1 alloc
                            #q_vec[j] += detJarray[i] * q_vec_i[j]  # -> 0 alloc
          end


        end # for i = 1 : integ_num

        
#=      for i = 1 : 8
            Q[ 1 + (elementmat[i,e]-1)*3 ] += q_vec[1 + (i-1)*3 ];
            Q[ 2 + (elementmat[i,e]-1)*3 ] += q_vec[2 + (i-1)*3 ];
            Q[ 3 + (elementmat[i,e]-1)*3 ] += q_vec[3 + (i-1)*3 ];
        end 
=#      
       
    end


#=    @inbounds for e = 1 : nElement
        @inbounds for i = 1 : 8
            Q[ 1 + (elementmat[i,e]-1)*3 ] += q_mat[1 + (i-1)*3, e];
            Q[ 2 + (elementmat[i,e]-1)*3 ] += q_mat[2 + (i-1)*3, e];
            Q[ 3 + (elementmat[i,e]-1)*3 ] += q_mat[3 + (i-1)*3, e];
        end
    end=#

            #integ_strain_ += d_integ_strain_  # -> same

            #return d_integ_stress_, d_integ_strain_, d_integ_yield_stress_, d_integ_plastic_strain_, d_integ_eq_plastic_strain_, Q
            #return d_integ_stress_, d_integ_strain_, d_integ_yield_stress_, d_integ_plastic_strain_, d_integ_eq_plastic_strain_  # -> OK
            #return d_integ_stress_
    return

end

#function cal_BVbar(Pusai_i, J_i, detJ_i)
#=function cal_BVbar(P2, detJ_i)

    BV_i = @MMatrix zeros(6,24);  # @MMatrix -> huge reduce alloc and memory
    #BVbar_i = zeros(6,24);

    #P2 = inv(J_i) * Pusai_i;
    #P2 = J_i \ Pusai_i;  # inv() is faster than \ in Julia in this case

    for i = 1:24
        BV_i[1,i] = P2[i]  # -> reduce alloc of N = reshape(P2,1,24)
        BV_i[2,i] = P2[i]
        BV_i[3,i] = P2[i] 
    end

#    N = reshape(P2,1,24);
#    BV_i[1,:] = N;
#    BV_i[2,:] = N;
#    BV_i[3,:] = N;

#=    for i = 1:24    # -> same as "BV_i[1,:] = N"
        BV_i[1,i] = N[i]
        BV_i[2,i] = N[i]
        BV_i[3,i] = N[i]
    end=#
    BV_i = BV_i / 3.0;
    #%BV_i = [N; N; N; zeros(1,24); zeros(1,24); zeros(1,24)] / 3;

    BVbar_i::MMatrix{6, 24, Float64, 144} = BV_i * detJ_i; # ::MMatrix{6, 24, Float64, 144} -> even

    #BVbar_i = @MMatrix zeros(6,24); #-> slow
    #BVbar_i .= BV_i * detJ_i;

    return BV_i, BVbar_i

end
=#


#function cal_Bfinal(Bfinal::MMatrix{6,24,Float64,144},
#                    Parray_i::MMatrix{3,8,Float64,24}, 
#                    BVbar::MMatrix{6,24,Float64,144})
function cal_Bfinal(Bfinal::MMatrix{6,24,Float64,144},
                    BVbar::MMatrix{6,24,Float64,144},
                    Pusai1::MMatrix{3,8,Float64,24}, 
                    e_position::MMatrix{3,8,Float64,24})

    J11 = J12 = J13 = 0.0
    J21 = J22 = J23 = 0.0
    J31 = J32 = J33 = 0.0

    @inbounds for i = 1 : 8
        J11 += Pusai1[1,i] * e_position[1,i]
        J12 += Pusai1[1,i] * e_position[2,i]
        J13 += Pusai1[1,i] * e_position[3,i]
        J21 += Pusai1[2,i] * e_position[1,i]
        J22 += Pusai1[2,i] * e_position[2,i]
        J23 += Pusai1[2,i] * e_position[3,i]
        J31 += Pusai1[3,i] * e_position[1,i]
        J32 += Pusai1[3,i] * e_position[2,i]
        J33 += Pusai1[3,i] * e_position[3,i]
    end

    v = ( J11*J22*J33 
        + J12*J23*J31 
        + J13*J21*J32 
        - J11*J23*J32 
        - J12*J21*J33 
        - J13*J22*J31 )
    detJi = v
    div_v::Float64 = 1.0/v

    iJ11 = ( J22*J33 - J23*J32 ) * div_v  # / v;
    iJ21 = ( J23*J31 - J21*J33 ) * div_v  # / v;
    iJ31 = ( J21*J32 - J22*J31 ) * div_v  # / v;

    iJ12 = ( J13*J32 - J12*J33 ) * div_v  # / v;
    iJ22 = ( J11*J33 - J13*J31 ) * div_v  # / v;
    iJ32 = ( J12*J31 - J11*J32 ) * div_v  # / v;

    iJ13 = ( J12*J23 - J13*J22 ) * div_v  # / v;
    iJ23 = ( J13*J21 - J11*J23 ) * div_v  # / v;
    iJ33 = ( J11*J22 - J12*J21 ) * div_v  # / v;


    @inbounds for i = 1 : 8

        #P2[1,i] = iJ11 * Pusai1[1,i] + iJ12 * Pusai1[2,i] + iJ13 * Pusai1[3,i] 
        #P2[2,i] = iJ21 * Pusai1[1,i] + iJ22 * Pusai1[2,i] + iJ23 * Pusai1[3,i] 
        #P2[3,i] = iJ31 * Pusai1[1,i] + iJ32 * Pusai1[2,i] + iJ33 * Pusai1[3,i] 

        Pix = iJ11 * Pusai1[1,i] + iJ12 * Pusai1[2,i] + iJ13 * Pusai1[3,i] 
        Piy = iJ21 * Pusai1[1,i] + iJ22 * Pusai1[2,i] + iJ23 * Pusai1[3,i] 
        Piz = iJ31 * Pusai1[1,i] + iJ32 * Pusai1[2,i] + iJ33 * Pusai1[3,i] 

        #Pix = Parray_i[1,i]
        #Piy = Parray_i[2,i] 
        #Piz = Parray_i[3,i]

        Bfinal[1,(i-1)*3+1] += Pix
        Bfinal[2,(i-1)*3+2] += Piy
        Bfinal[3,(i-1)*3+3] += Piz
        Bfinal[4,(i-1)*3+1] += Piy
        Bfinal[4,(i-1)*3+2] += Pix
        Bfinal[5,(i-1)*3+2] += Piz
        Bfinal[5,(i-1)*3+3] += Piy
        Bfinal[6,(i-1)*3+1] += Piz
        Bfinal[6,(i-1)*3+3] += Pix

        Bfinal[1,i*3-2] += -Pix  / 3.0  +  BVbar[1,i*3-2]
        Bfinal[1,i*3-1] += -Piy  / 3.0  +  BVbar[1,i*3-1]
        Bfinal[1,i*3]   += -Piz  / 3.0  +  BVbar[1,i*3]
        Bfinal[2,i*3-2] += -Pix  / 3.0  +  BVbar[2,i*3-2]
        Bfinal[2,i*3-1] += -Piy  / 3.0  +  BVbar[2,i*3-1]
        Bfinal[2,i*3]   += -Piz  / 3.0  +  BVbar[2,i*3]
        Bfinal[3,i*3-2] += -Pix  / 3.0  +  BVbar[3,i*3-2]
        Bfinal[3,i*3-1] += -Piy  / 3.0  +  BVbar[3,i*3-1]
        Bfinal[3,i*3]   += -Piz  / 3.0  +  BVbar[3,i*3]

#=      Bfinal[1,i*3-2] -= Pix  / 3.0
        Bfinal[1,i*3-1] -= Piy  / 3.0
        Bfinal[1,i*3]   -= Piz  / 3.0
        Bfinal[2,i*3-2] -= Pix  / 3.0
        Bfinal[2,i*3-1] -= Piy  / 3.0
        Bfinal[2,i*3]   -= Piz  / 3.0
        Bfinal[3,i*3-2] -= Pix  / 3.0
        Bfinal[3,i*3-1] -= Piy  / 3.0
        Bfinal[3,i*3]   -= Piz  / 3.0

        Bfinal[1,i*3-2] += BVbar[1,i*3-2]
        Bfinal[1,i*3-1] += BVbar[1,i*3-1]
        Bfinal[1,i*3]   += BVbar[1,i*3]
        Bfinal[2,i*3-2] += BVbar[2,i*3-2]
        Bfinal[2,i*3-1] += BVbar[2,i*3-1]
        Bfinal[2,i*3]   += BVbar[2,i*3]
        Bfinal[3,i*3-2] += BVbar[3,i*3-2]
        Bfinal[3,i*3-1] += BVbar[3,i*3-1]
        Bfinal[3,i*3]   += BVbar[3,i*3]
        =#

    end

    #Bfinal .+= BVbar

    return detJi

end

#=
function cal_P_hexa(Pusai1, e_position, 
                    Parray_i::MMatrix{3,8,Float64,24} ) 

    J11 = J12 = J13 = 0.0
    J21 = J22 = J23 = 0.0
    J31 = J32 = J33 = 0.0

    @inbounds for i = 1 : 8
        J11 += Pusai1[1,i] * e_position[1,i]
        J12 += Pusai1[1,i] * e_position[2,i]
        J13 += Pusai1[1,i] * e_position[3,i]
        J21 += Pusai1[2,i] * e_position[1,i]
        J22 += Pusai1[2,i] * e_position[2,i]
        J23 += Pusai1[2,i] * e_position[3,i]
        J31 += Pusai1[3,i] * e_position[1,i]
        J32 += Pusai1[3,i] * e_position[2,i]
        J33 += Pusai1[3,i] * e_position[3,i]
    end

    v = ( J11*J22*J33 
        + J12*J23*J31 
        + J13*J21*J32 
        - J11*J23*J32 
        - J12*J21*J33 
        - J13*J22*J31 )
    detJi = v
    div_v::Float64 = 1.0/v

    iJ11 = ( J22*J33 - J23*J32 ) * div_v  # / v;
    iJ21 = ( J23*J31 - J21*J33 ) * div_v  # / v;
    iJ31 = ( J21*J32 - J22*J31 ) * div_v  # / v;

    iJ12 = ( J13*J32 - J12*J33 ) * div_v  # / v;
    iJ22 = ( J11*J33 - J13*J31 ) * div_v  # / v;
    iJ32 = ( J12*J31 - J11*J32 ) * div_v  # / v;

    iJ13 = ( J12*J23 - J13*J22 ) * div_v  # / v;
    iJ23 = ( J13*J21 - J11*J23 ) * div_v  # / v;
    iJ33 = ( J11*J22 - J12*J21 ) * div_v  # / v;

    
    @inbounds for i = 1 : 8
        #P2[1,i] = iJ11 * Pusai1[1,i] + iJ12 * Pusai1[2,i] + iJ13 * Pusai1[3,i] 
        #P2[2,i] = iJ21 * Pusai1[1,i] + iJ22 * Pusai1[2,i] + iJ23 * Pusai1[3,i] 
        #P2[3,i] = iJ31 * Pusai1[1,i] + iJ32 * Pusai1[2,i] + iJ33 * Pusai1[3,i] 

        Pix = iJ11 * Pusai1[1,i] + iJ12 * Pusai1[2,i] + iJ13 * Pusai1[3,i] 
        Piy = iJ21 * Pusai1[1,i] + iJ22 * Pusai1[2,i] + iJ23 * Pusai1[3,i] 
        Piz = iJ31 * Pusai1[1,i] + iJ32 * Pusai1[2,i] + iJ33 * Pusai1[3,i] 

        Parray_i[1,i] = Pix
        Parray_i[2,i] = Piy
        Parray_i[3,i] = Piz
    end

    return detJi

end
=#


#function cal_B_BVbar_hexa(Pusai1, e_position, B, BV_i, BVbar_i)
#function cal_B_BVbar_hexa(Pusai1, e_position, B::MMatrix{6,24,Float64,144}, 
#                                              BV_i::MMatrix{6,24,Float64,144}, 
#                                              BVbar::MMatrix{6,24,Float64,144}) # ::MMatrix{6,24,Float64,144} -> fast
#=
function cal_B_BVbar_hexa(Pusai1, e_position, 
                            B::MMatrix{6,24,Float64,144}, 
                            BV_i::MMatrix{6,24,Float64,144},
                            Parray_i::MMatrix{3,8,Float64,24} ) 
    
    J11 = J12 = J13 = 0.0
    J21 = J22 = J23 = 0.0
    J31 = J32 = J33 = 0.0

    @inbounds for i = 1 : 8
        J11 += Pusai1[1,i] * e_position[1,i]
        J12 += Pusai1[1,i] * e_position[2,i]
        J13 += Pusai1[1,i] * e_position[3,i]
        J21 += Pusai1[2,i] * e_position[1,i]
        J22 += Pusai1[2,i] * e_position[2,i]
        J23 += Pusai1[2,i] * e_position[3,i]
        J31 += Pusai1[3,i] * e_position[1,i]
        J32 += Pusai1[3,i] * e_position[2,i]
        J33 += Pusai1[3,i] * e_position[3,i]
    end

    v = ( J11*J22*J33 
        + J12*J23*J31 
        + J13*J21*J32 
        - J11*J23*J32 
        - J12*J21*J33 
        - J13*J22*J31 )
    detJi = v
    div_v::Float64 = 1.0/v

    iJ11 = ( J22*J33 - J23*J32 ) * div_v  # / v;
    iJ21 = ( J23*J31 - J21*J33 ) * div_v  # / v;
    iJ31 = ( J21*J32 - J22*J31 ) * div_v  # / v;

    iJ12 = ( J13*J32 - J12*J33 ) * div_v  # / v;
    iJ22 = ( J11*J33 - J13*J31 ) * div_v  # / v;
    iJ32 = ( J12*J31 - J11*J32 ) * div_v  # / v;

    iJ13 = ( J12*J23 - J13*J22 ) * div_v  # / v;
    iJ23 = ( J13*J21 - J11*J23 ) * div_v  # / v;
    iJ33 = ( J11*J22 - J12*J21 ) * div_v  # / v;

        #P2 = zeros(3,8)
        #P2 .= 0.0
    #B .= 0
    #B = @MMatrix zeros(6,24)
        #div3::Float64 = 1.0/3.0

    @inbounds for i = 1 : 8
        #P2[1,i] = iJ11 * Pusai1[1,i] + iJ12 * Pusai1[2,i] + iJ13 * Pusai1[3,i] 
        #P2[2,i] = iJ21 * Pusai1[1,i] + iJ22 * Pusai1[2,i] + iJ23 * Pusai1[3,i] 
        #P2[3,i] = iJ31 * Pusai1[1,i] + iJ32 * Pusai1[2,i] + iJ33 * Pusai1[3,i] 

        Pix = iJ11 * Pusai1[1,i] + iJ12 * Pusai1[2,i] + iJ13 * Pusai1[3,i] 
        Piy = iJ21 * Pusai1[1,i] + iJ22 * Pusai1[2,i] + iJ23 * Pusai1[3,i] 
        Piz = iJ31 * Pusai1[1,i] + iJ32 * Pusai1[2,i] + iJ33 * Pusai1[3,i] 

        Parray_i[1,i] = Pix
        Parray_i[2,i] = Piy
        Parray_i[3,i] = Piz

        B[1,(i-1)*3+1] = Pix
        B[2,(i-1)*3+2] = Piy
        B[3,(i-1)*3+3] = Piz
        B[4,(i-1)*3+1] = Piy
        B[4,(i-1)*3+2] = Pix
        B[5,(i-1)*3+2] = Piz
        B[5,(i-1)*3+3] = Piy
        B[6,(i-1)*3+1] = Piz
        B[6,(i-1)*3+3] = Pix

        BV_i[1,i*3-2] = Pix  / 3.0
        BV_i[1,i*3-1] = Piy  / 3.0
        BV_i[1,i*3]   = Piz  / 3.0
        BV_i[2,i*3-2] = Pix  / 3.0
        BV_i[2,i*3-1] = Piy  / 3.0
        BV_i[2,i*3]   = Piz  / 3.0
        BV_i[3,i*3-2] = Pix  / 3.0
        BV_i[3,i*3-1] = Piy  / 3.0
        BV_i[3,i*3]   = Piz  / 3.0

    end

    #=B .= 0
    @inbounds for i = 1 : 8
      B[1,(i-1)*3+1] = P2[1,i]
      B[2,(i-1)*3+2] = P2[2,i]
      B[3,(i-1)*3+3] = P2[3,i]
      B[4,(i-1)*3+1] = P2[2,i]
      B[4,(i-1)*3+2] = P2[1,i]
      B[5,(i-1)*3+2] = P2[3,i]
      B[5,(i-1)*3+3] = P2[2,i]
      B[6,(i-1)*3+1] = P2[3,i]
      B[6,(i-1)*3+3] = P2[1,i]
    end=#

#=    @inbounds for i = 1:24
        BV_i[1,i] = P2[i] / 3.0  # -> reduce alloc of N = reshape(P2,1,24)
        BV_i[2,i] = P2[i] / 3.0
        BV_i[3,i] = P2[i] / 3.0
        # BV_i[4,i] = BV_i[5,i] = BV_i[6,i] = 0.0
    end =#

#=  @inbounds for i = 1:24
        #BVbar_i[1,i] = BV_i[1,i] * detJi
        #BVbar_i[2,i] = BV_i[2,i] * detJi
        #BVbar_i[3,i] = BV_i[3,i] * detJi
        BVbar[1,i] += BV_i[1,i] * detJi
        BVbar[2,i] += BV_i[2,i] * detJi
        BVbar[3,i] += BV_i[3,i] * detJi
    end  =#

    return detJi
    #return detJi, B
end
=#

function cal_BVbar_hexa(Pusai_mat, 
                        e_position::MMatrix{3,8,Float64,24}, 
                        BVbar::MMatrix{6,24,Float64,144}) # ::MMatrix{6,24,Float64,144} -> fast

    V = 0.0
    @inbounds for k = 1 : 8 #integ_num

        J11 = J12 = J13 = 0.0
        J21 = J22 = J23 = 0.0
        J31 = J32 = J33 = 0.0
        Pusai1 = Pusai_mat[k]

        @inbounds for i = 1 : 8
            J11 += Pusai1[1,i] * e_position[1,i]
            J12 += Pusai1[1,i] * e_position[2,i]
            J13 += Pusai1[1,i] * e_position[3,i]
            J21 += Pusai1[2,i] * e_position[1,i]
            J22 += Pusai1[2,i] * e_position[2,i]
            J23 += Pusai1[2,i] * e_position[3,i]
            J31 += Pusai1[3,i] * e_position[1,i]
            J32 += Pusai1[3,i] * e_position[2,i]
            J33 += Pusai1[3,i] * e_position[3,i]
        end

        v = ( J11*J22*J33 
        + J12*J23*J31 
        + J13*J21*J32 
        - J11*J23*J32 
        - J12*J21*J33 
        - J13*J22*J31 )
        detJi = v
        if detJi < 0
            detJi = abs(detJi)
            println("Warning:Element volume negative")
        end
        V += detJi

        div_v::Float64 = 1.0/detJi
        
        iJ11 = ( J22*J33 - J23*J32 ) * div_v  # / v;
        iJ21 = ( J23*J31 - J21*J33 ) * div_v  # / v;
        iJ31 = ( J21*J32 - J22*J31 ) * div_v  # / v;

        iJ12 = ( J13*J32 - J12*J33 ) * div_v  # / v;
        iJ22 = ( J11*J33 - J13*J31 ) * div_v  # / v;
        iJ32 = ( J12*J31 - J11*J32 ) * div_v  # / v;

        iJ13 = ( J12*J23 - J13*J22 ) * div_v  # / v;
        iJ23 = ( J13*J21 - J11*J23 ) * div_v  # / v;
        iJ33 = ( J11*J22 - J12*J21 ) * div_v  # / v;

        
        @inbounds for i = 1 : 8
            #P2[1,i] = iJ11 * Pusai1[1,i] + iJ12 * Pusai1[2,i] + iJ13 * Pusai1[3,i] 
            #P2[2,i] = iJ21 * Pusai1[1,i] + iJ22 * Pusai1[2,i] + iJ23 * Pusai1[3,i] 
            #P2[3,i] = iJ31 * Pusai1[1,i] + iJ32 * Pusai1[2,i] + iJ33 * Pusai1[3,i] 

            Pix = iJ11 * Pusai1[1,i] + iJ12 * Pusai1[2,i] + iJ13 * Pusai1[3,i] 
            Piy = iJ21 * Pusai1[1,i] + iJ22 * Pusai1[2,i] + iJ23 * Pusai1[3,i] 
            Piz = iJ31 * Pusai1[1,i] + iJ32 * Pusai1[2,i] + iJ33 * Pusai1[3,i] 
            
            BVbar[1,i*3-2] += Pix  / 3.0 * detJi
            BVbar[1,i*3-1] += Piy  / 3.0 * detJi
            BVbar[1,i*3]   += Piz  / 3.0 * detJi
            BVbar[2,i*3-2] += Pix  / 3.0 * detJi
            BVbar[2,i*3-1] += Piy  / 3.0 * detJi
            BVbar[2,i*3]   += Piz  / 3.0 * detJi
            BVbar[3,i*3-2] += Pix  / 3.0 * detJi
            BVbar[3,i*3-1] += Piy  / 3.0 * detJi
            BVbar[3,i*3]   += Piz  / 3.0 * detJi

        end

    end

    BVbar .= BVbar / V

    return V

end


#=
function cal_BVbar(P2, detJ_i, BV_i, BVbar_i)

    #BV_i = @MMatrix zeros(6,24);  # @MMatrix -> huge reduce alloc and memory

    @inbounds for i = 1:24
        BV_i[1,i] = P2[i] / 3.0  # -> reduce alloc of N = reshape(P2,1,24)
        BV_i[2,i] = P2[i] / 3.0
        BV_i[3,i] = P2[i] / 3.0
        # BV_i[4,i] = BV_i[5,i] = BV_i[6,i] = 0.0
    end

        #BV_i .= BV_i / 3.0;
                #BVbar_i::MMatrix{6, 24, Float64, 144} = BV_i * detJ_i; # ::MMatrix{6, 24, Float64, 144} -> even
    #BVbar_i .= BV_i * detJ_i
    @inbounds for i = 1:24
        BVbar_i[1,i] = BV_i[1,i] * detJ_i
        BVbar_i[2,i] = BV_i[2,i] * detJ_i
        BVbar_i[3,i] = BV_i[3,i] * detJ_i
    end

            #return BVbar_i
    return

end
=#


#=
function cal_B_hexa( Pusai1, e_position, B, P2 )
#function cal_B_hexa( Pusai1::Array{Float64,2}, e_position::MMatrix{3,8,Float64,24}, B::MMatrix{6,24,Float64,144}, P2::MMatrix{3,8,Float64,24} ) # no effect
    #function cal_B_hexa( Pusai1, e_position, B, P2, detJarray, index )

            #J = Pusai1 * e_position
        #J = Pusai1 * e_position'  # 0.8M/20k alloc
            #P2 .= 0.0
        #P2 .= my3inv(J) * Pusai1; # 1.6M/20k alloc
            #P2 = inv(J) * Pusai1;
            #P2 = J \ Pusai1;

    J11 = J12 = J13 = 0.0
    J21 = J22 = J23 = 0.0
    J31 = J32 = J33 = 0.0

    @inbounds for i = 1 : 8
        J11 += Pusai1[1,i] * e_position[1,i]
        J12 += Pusai1[1,i] * e_position[2,i]
        J13 += Pusai1[1,i] * e_position[3,i]
        J21 += Pusai1[2,i] * e_position[1,i]
        J22 += Pusai1[2,i] * e_position[2,i]
        J23 += Pusai1[2,i] * e_position[3,i]
        J31 += Pusai1[3,i] * e_position[1,i]
        J32 += Pusai1[3,i] * e_position[2,i]
        J33 += Pusai1[3,i] * e_position[3,i]
    end

    v = ( J11*J22*J33 
        + J12*J23*J31 
        + J13*J21*J32 
        - J11*J23*J32 
        - J12*J21*J33 
        - J13*J22*J31 )

    iJ11 = ( J22*J33 - J23*J32 ) / v;
    iJ21 = ( J23*J31 - J21*J33 ) / v;
    iJ31 = ( J21*J32 - J22*J31 ) / v;

    iJ12 = ( J13*J32 - J12*J33 ) / v;
    iJ22 = ( J11*J33 - J13*J31 ) / v;
    iJ32 = ( J12*J31 - J11*J32 ) / v;

    iJ13 = ( J12*J23 - J13*J22 ) / v;
    iJ23 = ( J13*J21 - J11*J23 ) / v;
    iJ33 = ( J11*J22 - J12*J21 ) / v;

    P2 .= 0.0
    @inbounds for i = 1 : 8
        P2[1,i] = iJ11 * Pusai1[1,i] + iJ12 * Pusai1[2,i] + iJ13 * Pusai1[3,i] 
        P2[2,i] = iJ21 * Pusai1[1,i] + iJ22 * Pusai1[2,i] + iJ23 * Pusai1[3,i] 
        P2[3,i] = iJ31 * Pusai1[1,i] + iJ32 * Pusai1[2,i] + iJ33 * Pusai1[3,i] 
    end

        #B = @MMatrix zeros(6,24);
    B .= 0
    @inbounds for i = 1 : 8
      B[1,(i-1)*3+1] = P2[1,i]
      B[2,(i-1)*3+2] = P2[2,i]
      B[3,(i-1)*3+3] = P2[3,i]
      B[4,(i-1)*3+1] = P2[2,i]
      B[4,(i-1)*3+2] = P2[1,i]
      B[5,(i-1)*3+2] = P2[3,i]
      B[5,(i-1)*3+3] = P2[2,i]
      B[6,(i-1)*3+1] = P2[3,i]
      B[6,(i-1)*3+3] = P2[1,i]
    end

        #detJi = 0.0
        #detJi = my3det(J)
    detJi = v
    #detJarray[index] = v

        #return B, J, P2
        #return B, P2, detJi
    return detJi

end
=#

function cal_Pusai_hexa(integ_num)

    #Pusai_mat = zeros(3,8,integ_num); #@MArray -> huge increase alloc and memory
    Pusai_mat = [ @MMatrix zeros(3,8) for i=1:8]  # -> huge reduce alloc and memory compared to  zeros(3,8,integ_num)
    
    delta_mat = [ -1.0  -1.0  -1.0
                   1.0  -1.0  -1.0
                   1.0   1.0  -1.0
                  -1.0   1.0  -1.0
                  -1.0  -1.0   1.0
                   1.0  -1.0   1.0
                   1.0   1.0   1.0
                  -1.0   1.0   1.0 ];

     # integral point
     gc = [0 0 0]; 
     if integ_num == 8
         g = 1.0 / sqrt(3.0);
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
     Pusai1 = @MMatrix zeros(3,8); #@MMatrix -> even 
     gzai = gc[k,1];
     eta = gc[k,2];
     tueta = gc[k,3];

      for i = 1 : 8
        Pusai1[1,i] = 1.0/8.0*delta_mat[i,1]*(1.0+eta* delta_mat[i,2])*(1.0+tueta* delta_mat[i,3]);
        Pusai1[2,i] = 1.0/8.0*delta_mat[i,2]*(1.0+gzai* delta_mat[i,1])*(1.0+tueta* delta_mat[i,3]);
        Pusai1[3,i] = 1.0/8.0*delta_mat[i,3]*(1.0+gzai* delta_mat[i,1])*(1.0+eta* delta_mat[i,2]);
      end

      #Pusai_mat[:,:,k] = Pusai1;
      Pusai_mat[k] .= Pusai1;

    end

    return Pusai_mat

end


function get_element_face(MODEL, i)

    part_id = MODEL.INSTANCE[i].part_id;
    cdmat = MODEL.PART[ part_id ].coordmat;
    nE = MODEL.INSTANCE[i].nElement;
    faces = zeros(Int, nE*6, 4);
    faces_eleid = zeros(Int, nE*6);
    sorted_faces = zeros(Int, nE*6, 4);

    c = 1;
    for j = 1 : nE
        #elem = MODEL.PART[ part_id ].elementmat[j,:];
        elem = MODEL.PART[ part_id ].elementmat[:,j];
        faces[6*(c-1)+1, :] = elem[1:4];
        faces[6*(c-1)+2, :] = elem[5:8];
        faces[6*(c-1)+3, :] = [elem[1], elem[2], elem[6], elem[5]];
        faces[6*(c-1)+4, :] = [elem[2], elem[3], elem[7], elem[6]];
        faces[6*(c-1)+5, :] = [elem[3], elem[4], elem[8], elem[7]];
        faces[6*(c-1)+6, :] = [elem[4], elem[1], elem[5], elem[8]];
        faces_eleid[6*(c-1).+(1:6)] .= j;
        #ctr = sum( cdmat[elem,:] ) / 8;
        #ctr = sum!(zeros(1,3), cdmat[elem,:] ) / 8
        ctr = sum!(zeros(3,1), cdmat[:,elem] ) / 8
        #println("ctr:",ctr)
        for k = 1 : 6
            index = 6*(c-1)+k;
            #v1 = cdmat[ faces[index, 2], :] - cdmat[ faces[index, 1], :];
            #v2 = cdmat[ faces[index, 4], :] - cdmat[ faces[index, 1], :];
            v1 = cdmat[ :, faces[index, 2]] - cdmat[ :, faces[index, 1]];
            v2 = cdmat[ :, faces[index, 4]] - cdmat[ :, faces[index, 1]];
            nv = my3cross(v1,v2);
            #vc = ctr' - cdmat[ faces[index, 1], :];
            vc = ctr - cdmat[:, faces[index, 1]];
            if dot(nv,vc) > 0.
                #%'order modified'
                faces[index, :] = [faces[index, 1],  faces[index, 4], faces[index, 3], faces[index, 2]];
            end
        end
        c = c + 1;
    end

    for j = 1 : nE*6
        sorted_faces[j,:] = sort(faces[j,:]);
    end

    return faces, faces_eleid, sorted_faces
end


#function get_surface_triangle(MODEL, i, array_element)
function get_surface_triangle(INSTANCE_i, array_element, contact_element)

    nE = length(array_element); 
    surfaces = zeros(Int, nE*6, 4);
    surfaces_eleid = zeros(Int, nE*6);
    sorted_surfaces = zeros(Int, nE*6, 4);
    c = 1;
    for j = array_element
        for k = 1 : 6
                        #surfaces[6*(c-1)+k, :] = MODEL.INSTANCE[i].surfaces[6*(j-1)+k,:];
                        #sorted_surfaces[6*(c-1)+k, :] = MODEL.INSTANCE[i].sorted_surfaces[6*(j-1)+k,:];
                #surfaces[6*(c-1)+k, :] = INSTANCE_i.surfaces[6*(j-1)+k,:];  #->alloc
                #sorted_surfaces[6*(c-1)+k, :] = INSTANCE_i.sorted_surfaces[6*(j-1)+k,:]; #->alloc
            surfaces[6*(c-1)+k, 1] = INSTANCE_i.surfaces[6*(j-1)+k, 1]
            surfaces[6*(c-1)+k, 2] = INSTANCE_i.surfaces[6*(j-1)+k, 2]
            surfaces[6*(c-1)+k, 3] = INSTANCE_i.surfaces[6*(j-1)+k, 3]
            surfaces[6*(c-1)+k, 4] = INSTANCE_i.surfaces[6*(j-1)+k, 4]
            sorted_surfaces[6*(c-1)+k, 1] = INSTANCE_i.sorted_surfaces[6*(j-1)+k, 1]
            sorted_surfaces[6*(c-1)+k, 2] = INSTANCE_i.sorted_surfaces[6*(j-1)+k, 2]
            sorted_surfaces[6*(c-1)+k, 3] = INSTANCE_i.sorted_surfaces[6*(j-1)+k, 3]
            sorted_surfaces[6*(c-1)+k, 4] = INSTANCE_i.sorted_surfaces[6*(j-1)+k, 4]
        end
        
        surfaces_eleid[6*(c-1).+(1:6)] .= j;
            #surfaces_eleid[6*(c-1)+1] = j   # -> same as above
            #surfaces_eleid[6*(c-1)+2] = j
            #surfaces_eleid[6*(c-1)+3] = j
            #surfaces_eleid[6*(c-1)+4] = j
            #surfaces_eleid[6*(c-1)+5] = j
            #surfaces_eleid[6*(c-1)+6] = j

        c = c + 1;
    end


    #c_surfaces = Vector{Vector{Int}}();
    #c_surfaces_eleid = Int[];
    c_surfaces = zeros(Int, nE*6, 4);
    c_surfaces_eleid = zeros(Int, nE*6);
    dp_id = Int[];
    sj = zeros(Int,4)
    sk = zeros(Int,4)
    c = 0
    #@time 
    for j = 1 : nE*6-1
        u_flag = 1;
        #if length(findall(dp_id.==j)) > 0
        if !isnothing(findfirst(dp_id.==j)) 
            u_flag = 0;
            continue;
        end

        #sj = sorted_surfaces[j, :]
        sj[1] = sorted_surfaces[j, 1]
        sj[2] = sorted_surfaces[j, 2]
        sj[3] = sorted_surfaces[j, 3]
        sj[4] = sorted_surfaces[j, 4]
        
        for k = j+1 : nE*6  #%1 : nE*6

            #sk = sorted_surfaces[k, :]
            sk[1] = sorted_surfaces[k, 1]
            sk[2] = sorted_surfaces[k, 2]
            sk[3] = sorted_surfaces[k, 3]
            sk[4] = sorted_surfaces[k, 4]

            if sj[1] == sk[1] && sj[2] == sk[2] && sj[3] == sk[3] && sj[4] == sk[4]
                #sj == sk # -> slow

                u_flag = 0;
                #dp_id = [dp_id k];
                append!(dp_id, k)
                break
            end
        end

        if u_flag == 1
            #push!(c_surfaces, surfaces[j, :])
            #append!(c_surfaces_eleid, surfaces_eleid[j])
            c += 1
            #c_surfaces[c,:] = surfaces[j, :]
            c_surfaces[c,1] = surfaces[j,1]
            c_surfaces[c,2] = surfaces[j,2]
            c_surfaces[c,3] = surfaces[j,3]
            c_surfaces[c,4] = surfaces[j,4]
            c_surfaces_eleid[c] = surfaces_eleid[j]
        end

    end
    #length(c_surfaces);

    c_surfaces = c_surfaces[1:c,:]
    c_surfaces_eleid = c_surfaces_eleid[1:c]
    #println("c:", c)
    #println("c_surfaces:", c_surfaces)
    #println("c_surfaces_eleid:", c_surfaces_eleid)
    
    #%--- Pick up only contact element
    if INSTANCE_i.nElement != length(contact_element)
        #c_surfaces_temp = Vector{Vector{Int}}()
        #c_surfaces_eleid_temp = Int[];
        c_surfaces_temp = zeros(Int, c, 4)
        c_surfaces_eleid_temp = zeros(Int, c)
        c = 0
        for j = 1 : length(c_surfaces_eleid)
            #if length( findall(contact_element .== c_surfaces_eleid[j]) ) > 0
            if !isnothing( findfirst( contact_element .== c_surfaces_eleid[j] ) ) 
                #push!(c_surfaces_temp, c_surfaces[j])
                #append!(c_surfaces_eleid_temp, c_surfaces_eleid[j])
                c += 1
                #c_surfaces_temp[c,:] = c_surfaces[j, :]
                c_surfaces_temp[c,1] = c_surfaces[j,1]
                c_surfaces_temp[c,2] = c_surfaces[j,2]
                c_surfaces_temp[c,3] = c_surfaces[j,3]
                c_surfaces_temp[c,4] = c_surfaces[j,4]
                c_surfaces_eleid_temp[c] = c_surfaces_eleid[j]
            end
        end

        #c_surfaces = c_surfaces_temp;
        #c_surfaces_eleid = c_surfaces_eleid_temp;
        c_surfaces = c_surfaces_temp[1:c,:]
        c_surfaces_eleid = c_surfaces_eleid_temp[1:c]
    end


    if length(c_surfaces) == 0
        c_triangles = zeros(Int,0,0);
        c_triangles_eleid = zeros(Int,0);
        c_nodes = zeros(Int,0);
        return
    end

    c_triangles = zeros(Int, size(c_surfaces,1)*2,3);
    c_triangles_eleid = zeros(Int, size(c_surfaces,1)*2);
    for j = 1 : size(c_surfaces,1)
                #c_triangles[j*2-1,:] = [c_surfaces[j][1], c_surfaces[j][2], c_surfaces[j][3]];
                #c_triangles[j*2,:]   = [c_surfaces[j][3], c_surfaces[j][4], c_surfaces[j][1]];
            #c_triangles[j*2-1,1] = c_surfaces[j][1]
            #c_triangles[j*2-1,2] = c_surfaces[j][2]
            #c_triangles[j*2-1,3] = c_surfaces[j][3]
            #c_triangles[j*2,1] = c_surfaces[j][3]
            #c_triangles[j*2,2] = c_surfaces[j][4]
            #c_triangles[j*2,3] = c_surfaces[j][1]
        c_triangles[j*2-1,1] = c_surfaces[j,1]
        c_triangles[j*2-1,2] = c_surfaces[j,2]
        c_triangles[j*2-1,3] = c_surfaces[j,3]
        c_triangles[j*2,1] = c_surfaces[j,3]
        c_triangles[j*2,2] = c_surfaces[j,4]
        c_triangles[j*2,3] = c_surfaces[j,1]

        c_triangles_eleid[j*2-1] = c_surfaces_eleid[j];
        c_triangles_eleid[j*2] = c_surfaces_eleid[j];
    end
    
                #c_nodes = reshape(c_surfaces, size(c_surfaces,1)*4);
        #c_nodes = reduce(vcat, c_surfaces);
    c_nodes = zeros(Int, size(c_triangles,1)*3)
    for i = 1 : size(c_triangles,1)*3
        c_nodes[i] = c_triangles[i]
    end

        #c_nodes = sort(unique(c_nodes));
    sort!(unique!(c_nodes))
    #%length(c_nodes);
    
    return c_triangles, c_triangles_eleid, c_nodes

end


function add_surface_triangle(INSTANCE_i, ele_id::Int)

    # find surfaces which are match to the deleted element

    nE = INSTANCE_i.nElement
    surfaces = zeros(Int, 6, 4);
    sorted_surfaces = zeros(Int, 6, 4);
    for k = 1 : 6
        surfaces[k, 1] = INSTANCE_i.surfaces[6*(ele_id-1)+k, 1]
        surfaces[k, 2] = INSTANCE_i.surfaces[6*(ele_id-1)+k, 2]
        surfaces[k, 3] = INSTANCE_i.surfaces[6*(ele_id-1)+k, 3]
        surfaces[k, 4] = INSTANCE_i.surfaces[6*(ele_id-1)+k, 4]
        sorted_surfaces[k, 1] = INSTANCE_i.sorted_surfaces[6*(ele_id-1)+k, 1]
        sorted_surfaces[k, 2] = INSTANCE_i.sorted_surfaces[6*(ele_id-1)+k, 2]
        sorted_surfaces[k, 3] = INSTANCE_i.sorted_surfaces[6*(ele_id-1)+k, 3]
        sorted_surfaces[k, 4] = INSTANCE_i.sorted_surfaces[6*(ele_id-1)+k, 4]
    end

    
    add_c_surfaces = Vector{Vector{Int}}();
    add_c_surfaces_eleid = Int[];
    sj = zeros(Int,4)
    sk = zeros(Int,4)

    for j = 1 : 6
        
        #sj = sorted_surfaces[j, :]
        sj[1] = sorted_surfaces[j, 1]
        sj[2] = sorted_surfaces[j, 2]
        sj[3] = sorted_surfaces[j, 3]
        sj[4] = sorted_surfaces[j, 4]
        
        for k = 1 : nE*6  
            if INSTANCE_i.surfaces_eleid[k] == ele_id
                continue
            end

            #sk = sorted_surfaces[k, :]
            sk[1] = INSTANCE_i.sorted_surfaces[k, 1]
            sk[2] = INSTANCE_i.sorted_surfaces[k, 2]
            sk[3] = INSTANCE_i.sorted_surfaces[k, 3]
            sk[4] = INSTANCE_i.sorted_surfaces[k, 4]

            if sj[1] == sk[1] && sj[2] == sk[2] && sj[3] == sk[3] && sj[4] == sk[4]
                #sj == sk # -> slow

                push!(add_c_surfaces, INSTANCE_i.surfaces[k, :])
                append!(add_c_surfaces_eleid, INSTANCE_i.surfaces_eleid[k])
                break
            end
        end

    end    


    add_c_triangles = zeros(Int, size(add_c_surfaces,1)*2,3);
    add_c_triangles_eleid = zeros(Int, size(add_c_surfaces,1)*2);
    for j = 1 : size(add_c_surfaces,1)
        add_c_triangles[j*2-1,1] = add_c_surfaces[j][1]
        add_c_triangles[j*2-1,2] = add_c_surfaces[j][2]
        add_c_triangles[j*2-1,3] = add_c_surfaces[j][3]
        add_c_triangles[j*2,1] = add_c_surfaces[j][3]
        add_c_triangles[j*2,2] = add_c_surfaces[j][4]
        add_c_triangles[j*2,3] = add_c_surfaces[j][1]

        add_c_triangles_eleid[j*2-1] = add_c_surfaces_eleid[j];
        add_c_triangles_eleid[j*2] = add_c_surfaces_eleid[j];
    end
    
    add_c_nodes = zeros(Int, size(add_c_triangles,1)*3)
    for i = 1 : size(add_c_triangles,1)*3
        add_c_nodes[i] = add_c_triangles[i]
    end

    sort!(unique!(add_c_nodes))

    return add_c_triangles, add_c_triangles_eleid, add_c_nodes

end


function cal_contact_force(c_force3, CT, instance_pair, cp_index, # CV, 
                           position::Array{Float64,2}, velo, diag_M, elementMinSize, elementMaxSize, d_max,  
                           element_flag, elementmat::Array{Int,2}, bug_report, time_)
      
    nNode = size(position,2);

    d_lim::Float64 = elementMinSize * 0.3;
    myu = 0.25 * 1.0
    kc_o = 1.0;
    kc_s = 1.0; # self-contact
    Cr_o = 0.0;
    Cr_s = 0.0;

    div3::Float64 = 1.0/3.0


    for c = 1 : length(cp_index)

        cc = cp_index[c]
        
        #% j -> triangle,   i -> point
        i_instance = instance_pair[c][1]
        j_instance = instance_pair[c][2]
        #println("i_instance:", i_instance, " j_instance:", j_instance)
     
        c_nodes_i::Array{Int,1} = CT[c].c_nodes_i    # -> 0 alloc
        c_nodes_j::Array{Int,1} = CT[c].c_nodes_j
        c_triangles::Array{Int,2} = CT[c].c_triangles
        c_triangles_eleid::Array{Int,1} = CT[c].c_triangles_eleid
        young = CT[c].young
        nn_i = length(c_nodes_i)
        nn_j = length(c_nodes_j)
        #println("c_nodes_i:", length(c_nodes_i), ", c_triangles:", size(c_triangles,1) )


        #%--- contact range ---%    
        min_ix = minimum(@view(position[1,c_nodes_i]))
        min_iy = minimum(@view(position[2,c_nodes_i]))
        min_iz = minimum(@view(position[3,c_nodes_i]))
        min_jx = minimum(@view(position[1,c_nodes_j]))
        min_jy = minimum(@view(position[2,c_nodes_j]))
        min_jz = minimum(@view(position[3,c_nodes_j]))

        max_ix = maximum(@view(position[1,c_nodes_i]))
        max_iy = maximum(@view(position[2,c_nodes_i]))
        max_iz = maximum(@view(position[3,c_nodes_i]))
        max_jx = maximum(@view(position[1,c_nodes_j]))
        max_jy = maximum(@view(position[2,c_nodes_j]))
        max_jz = maximum(@view(position[3,c_nodes_j]))

        range_min_x = max(min_ix, min_jx)
        range_min_y = max(min_iy, min_jy)
        range_min_z = max(min_iz, min_jz)
        range_max_x = min(max_ix, max_jx)
        range_max_y = min(max_iy, max_jy)
        range_max_z = min(max_iz, max_jz)
    
        if range_min_x > range_max_x || range_min_y > range_max_y || range_min_z > range_max_z
            continue
        end


        all_range_min_x = min(min_ix, min_jx)
        all_range_min_y = min(min_iy, min_jy)
        all_range_min_z = min(min_iz, min_jz)
        all_range_max_x = max(max_ix, max_jx)
        all_range_max_y = max(max_iy, max_jy)
        all_range_max_z = max(max_iz, max_jz)
                
        RLx = all_range_max_x - all_range_min_x
        RLy = all_range_max_y - all_range_min_y
        RLz = all_range_max_z - all_range_min_z

        #nn_i = length(c_nodes_i)
        #nn_j = length(c_nodes_j)

        node_map_ix = zeros(Int, nn_i)  
        node_map_iy = zeros(Int, nn_i)  
        node_map_iz = zeros(Int, nn_i)  
        node_map_jx = zeros(Int, nn_j)  
        node_map_jy = zeros(Int, nn_j)  
        node_map_jz = zeros(Int, nn_j)  

        ddiv = elementMaxSize * 1.1
        if i_instance == j_instance
            ddiv = elementMaxSize * 0.6
        end

        
        @inbounds for i::Int = 1 : nn_i 
            pos_i_x = position[1, c_nodes_i[i]]
            map_i_x::Int = ceil( (pos_i_x - all_range_min_x) / ddiv )
            node_map_ix[i] = map_i_x

            pos_i_y = position[2, c_nodes_i[i]]
            map_i_y::Int = ceil( (pos_i_y - all_range_min_y) / ddiv )
            node_map_iy[i] = map_i_y

            pos_i_z = position[3, c_nodes_i[i]]
            map_i_z::Int = ceil( (pos_i_z - all_range_min_z) / ddiv )
            node_map_iz[i] = map_i_z
        end
        
        @inbounds for i::Int = 1 : nn_j
            pos_j_x = position[1, c_nodes_j[i]]
            map_j_x::Int = ceil( (pos_j_x - all_range_min_x) / ddiv )
            node_map_jx[i] = map_j_x

            pos_j_y = position[2, c_nodes_j[i]]
            map_j_y::Int = ceil( (pos_j_y - all_range_min_y) / ddiv )
            node_map_jy[i] = map_j_y

            pos_j_z = position[3, c_nodes_j[i]]
            map_j_z::Int = ceil( (pos_j_z - all_range_min_z) / ddiv )
            node_map_jz[i] = map_j_z
        end
        

        #array_randj = sortperm( rand(size(c_triangles,1)) )  # different result in self-contact
        #array_randj = shuffle(1 : size(c_triangles,1))  # different result in self-contact
        #println("sum(array_randj):", sum(array_randj))

        @inbounds @floop for j = 1 : size(c_triangles,1)
        #@inbounds @floop for j = array_randj

            eleid_ = c_triangles_eleid[j]
            if element_flag[ eleid_ ] == 0
                continue
            end

            kc = kc_o
            Cr = Cr_o
            if i_instance == j_instance
                kc = kc_s
                Cr = Cr_s
            end

            #index_th = 1 
            index_th = Threads.threadid()
            #l = ReentrantLock()  

            j0 = c_triangles[j,1]
            j1 = c_triangles[j,2]
            j2 = c_triangles[j,3]

            q0x = position[1,j0]
            q0y = position[2,j0]
            q0z = position[3,j0]
            q1x = position[1,j1]
            q1y = position[2,j1]
            q1z = position[3,j1]
            q2x = position[1,j2]
            q2y = position[2,j2]
            q2z = position[3,j2]

            if q0x < range_min_x && q1x < range_min_x && q2x < range_min_x
                continue
            end
            if q0y < range_min_y && q1y < range_min_y && q2y < range_min_y
                continue
            end
            if q0z < range_min_z && q1z < range_min_z && q2z < range_min_z
                continue
            end

            if q0x > range_max_x && q1x > range_max_x && q2x > range_max_x
                continue
            end
            if q0y > range_max_y && q1y > range_max_y && q2y > range_max_y
                continue
            end
            if q0z > range_max_z && q1z > range_max_z && q2z > range_max_z
                continue
            end

            cx = (q0x+q1x+q2x) /3.0  #* div3 # /3.0  #
            cy = (q0y+q1y+q2y) /3.0  #* div3 # /3.0
            cz = (q0z+q1z+q2z) /3.0  #* div3 # /3.0
            R0 = my3norm(q0x-cx, q0y-cy, q0z-cz)
            R1 = my3norm(q1x-cx, q1y-cy, q1z-cz)
            R2 = my3norm(q2x-cx, q2y-cy, q2z-cz)
            Rmax = max(max(R0,R1),R2)


            v1x = q1x - q0x
            v1y = q1y - q0y
            v1z = q1z - q0z

            v2x = q2x - q0x
            v2y = q2y - q0y
            v2z = q2z - q0z

            L1 = my3norm(v1x,v1y,v1z);
            L2 = my3norm(v2x,v2y,v2z);
            Lmax = max(L1,L2);
            nx, ny, nz = my3crossNNz(v1x,v1y,v1z, v2x,v2y,v2z)
            d12 = v1x*v2x+v1y*v2y+v1z*v2z
            S = 0.5 * sqrt( L1*L1 * L2*L2 - d12*d12 );

                    #--- custom
                    # if abs(n[3]) > 0.99
                    #     continue
                    # end

            A11 = v1x #v1[1] #->OK
            A21 = v1y #v1[2]
            A31 = v1z #v1[3]
            A12 = v2x #v2[1]
            A22 = v2y #v2[2]
            A32 = v2z #v2[3] #->OK
            A13 = -nx #-n[1]
            A23 = -ny #-n[2]
            A33 = -nz #-n[3]

            map_j0_x = 1
            map_j0_y = 1
            map_j0_z = 1
            @inbounds for i = 1 : nn_j
                if c_nodes_j[i] == j0
                    map_j0_x = node_map_jx[i]
                    map_j0_y = node_map_jy[i]
                    map_j0_z = node_map_jz[i]
                    break
                end
            end 

            e1_ = elementmat[1, eleid_]
            e2_ = elementmat[2, eleid_]
            e3_ = elementmat[3, eleid_]
            e4_ = elementmat[4, eleid_]
            e5_ = elementmat[5, eleid_]
            e6_ = elementmat[6, eleid_]
            e7_ = elementmat[7, eleid_]
            e8_ = elementmat[8, eleid_]


            #@inbounds for  i = c_nodes_i
            @inbounds for  k = 1 : nn_i
            
                if (abs(map_j0_x - node_map_ix[k]) > 1  ||
                    abs(map_j0_y - node_map_iy[k]) > 1  ||
                    abs(map_j0_z - node_map_iz[k]) > 1  )
                    
                    continue
                end 

                i = c_nodes_i[k]

                if i_instance == j_instance
                    if ( i == e1_ ||
                        i == e2_ ||
                        i == e3_ ||
                        i == e4_ ||
                        i == e5_ ||
                        i == e6_ ||
                        i == e7_ ||
                        i == e8_  )
                        continue
                    end #->no cost
                end

                px = position[1,i]
                py = position[2,i]
                pz = position[3,i]

                #print(bug_report, @sprintf("time=%.3e, i=%d, p=%.16e, %.16e, %.16e\n", 
                #                            time_, i, px, py, pz) ) 
      
                if px < range_min_x || py < range_min_y || pz < range_min_z
                    continue
                end
                if px > range_max_x || py > range_max_y || pz > range_max_z
                    continue
                end

                

                dpc = my3norm(px-cx,py-cy,pz-cz)
                if dpc >= Rmax
                    continue
                end

                bx = px - q0x
                by = py - q0y
                bz = pz - q0z

                #if my3norm(bx,by,bz) > Lmax
                #    continue
                #end

                        # dA = abs(det(A));
                        # if dA < 1E-10
                        #     dA
                        #     continue;
                        # end

                        # rA = abs(rcond(A));
                        # if rA < 1E-3 %|| isnan(rA)
                        #     rA
                        #     continue;
                        # end

                x1, x2, d = my3SolveAb(A11,A21,A31, A12,A22,A32, A13,A23,A33,
                                       bx,by,bz)

                if (0.0 <= x1 && 0.0 <= x2 && x1 + x2 <= 1.0 &&   # -> 10 alloc -> 0 alloc
                    d > 0.0 && d <= d_lim )

                    #=
                    if d - d_node_pre[i] > d_max
                        d = d_node_pre[i] + d_max;
                    end =#

                    vx = velo[i*3-2] - velo[j0*3-2]
                    vy = velo[i*3-1] - velo[j0*3-1]
                    vz = velo[i*3] - velo[j0*3]
                    mag_v = my3norm(vx,vy,vz)

                    vex = vey = vez = 0.0
                    
                    if mag_v > 0.0     # 1.0E-10  #0.0
                        vex = vx / mag_v
                        vey = vy / mag_v
                        vez = vz / mag_v  # -> 0 alloc
                            #cv = abs( ve[1]*n[1]+ve[2]*n[2]+ve[3]*n[3] )
                    end
                        #println("ve:", ve[1],",",ve[2],",",ve[3]," cv=", cv)
                        
                    k = young * S / Lmax  * kc  # * 2.0 ~ 10.0  # -> for hard contact
                                # cc = 2 * sqrt( diag_M(i) * k ) * 0;
                                #f = k * x[3] * n; #  - cc * (v*n.') * n;
                    F = k * d
                                #F = k * d * cv  
                    fx::Float64 = F * nx
                    fy::Float64 = F * ny
                    fz::Float64 = F * nz

                    #print(bug_report, @sprintf("time=%.3e, i=%d, cv=%.6e, vi=%.6e, %.6e, %.6e\n", 
                    #                            time_, i, cv, v[1], v[2], v[3]) )
                    #print(bug_report, @sprintf("time=%.3e, i=%d, vi=%.16e, %.16e, %.16e\n", 
                    #                            time_, i, v[1], v[2], v[3]) )
                    #print(bug_report, @sprintf("time=%.3e, i=%d, pi=%.16e, %.16e, %.16e\n", 
                    #                            time_, i, p[1], p[2], p[3]) )

                    # damping
                    C = 2 * sqrt( diag_M[i] * k ) * Cr;
                        #f = k * x[3] * n; #  - cc * (v*n.') * n;
                    fc_x = -C * vx
                    fc_y = -C * vy
                    fc_z = -C * vz


                    #%--- Friction
                            #vs .= ve .- ( ve[1]*n[1]+ve[2]*n[2]+ve[3]*n[3] ) * n  # much faster?  #->3 alloc
                        #dot_ve_n = ve[1]*n[1]+ve[2]*n[2]+ve[3]*n[3]
                        #vs[1] = ve[1] - dot_ve_n * n[1];  # become slow?
                        #vs[2] = ve[2] - dot_ve_n * n[2];
                        #vs[3] = ve[3] - dot_ve_n * n[3];  # -> 0 alloc
                    dot_ve_n = vex*nx + vey*ny + vez*nz
                    vsx = vex - dot_ve_n * nx  # become slow?
                    vsy = vey - dot_ve_n * ny
                    vsz = vez - dot_ve_n * nz  # -> 0 alloc
                   
                            #println("vs:", vs[1],",",vs[2],",",vs[3])
                            #println("fric:", fric[1], ",", fric[2], ",", fric[3])
                    fric_x::Float64 = -myu * F * vsx
                    fric_y::Float64 = -myu * F * vsy
                    fric_z::Float64 = -myu * F * vsz  # -> 0 alloc
                    fx += fric_x + fc_x
                    fy += fric_y + fc_y
                    fz += fric_z + fc_z
                
                    #print(bug_report, @sprintf("time=%.3e, i=%d, n=%.16e, %.16e, %.16e\n", 
                    #                            time_, i, n[1], n[2], n[3]) ) 
                    #print(bug_report, @sprintf("time=%.3e, i=%d, ve=%.16e, %.16e, %.16e\n", 
                    #                            time_, i, ve[1], ve[2], ve[3]) ) 
                    #print(bug_report, @sprintf("time=%.3e, i=%d, vs=%.16e, %.16e, %.16e\n", 
                    #                            time_, i, vs[1], vs[2], vs[3]) ) 

                    #println("f:", f[1], ", ", f[2], ", ", f[3])
                    #print(bug_report, @sprintf("time=%.3e, i=%d, k=%.16e, d=%.16e, S=%.16e\n", 
                    #                            time_, i, k, d, S) )  # -> all same 
                    #print(bug_report, @sprintf("time=%.3e, i=%d, f=%.16e, %.16e, %.16e\n", 
                    #                            time_, i, f[1], f[2], f[3]) ) #->diff
                    #print(bug_report, @sprintf("time=%.3e, i=%d, f=%.16e, %.16e, %.16e\n", 
                    #                            time_, i, fx, fy, fz) ) #->diff


                    #=  c_force3[1+(i-1)*3] +=  fx
                        c_force3[2+(i-1)*3] +=  fy
                        c_force3[3+(i-1)*3] +=  fz

                        c_force3[1+(j0-1)*3] +=  -fx  / 3.0
                        c_force3[2+(j0-1)*3] +=  -fy  / 3.0
                        c_force3[3+(j0-1)*3] +=  -fz  / 3.0

                        c_force3[1+(j1-1)*3] +=  -fx  / 3.0
                        c_force3[2+(j1-1)*3] +=  -fy  / 3.0
                        c_force3[3+(j1-1)*3] +=  -fz  / 3.0

                        c_force3[1+(j2-1)*3] +=  -fx  / 3.0
                        c_force3[2+(j2-1)*3] +=  -fy  / 3.0
                        c_force3[3+(j2-1)*3] +=  -fz  / 3.0
                    =#

                    c_force3[1+(i-1)*3, index_th] +=  fx
                    c_force3[2+(i-1)*3, index_th] +=  fy
                    c_force3[3+(i-1)*3, index_th] +=  fz

                    c_force3[1+(j0-1)*3, index_th] +=  -fx  / 3.0
                    c_force3[2+(j0-1)*3, index_th] +=  -fy  / 3.0
                    c_force3[3+(j0-1)*3, index_th] +=  -fz  / 3.0

                    c_force3[1+(j1-1)*3, index_th] +=  -fx  / 3.0
                    c_force3[2+(j1-1)*3, index_th] +=  -fy  / 3.0
                    c_force3[3+(j1-1)*3, index_th] +=  -fz  / 3.0

                    c_force3[1+(j2-1)*3, index_th] +=  -fx  / 3.0
                    c_force3[2+(j2-1)*3, index_th] +=  -fy  / 3.0
                    c_force3[3+(j2-1)*3, index_th] +=  -fz  / 3.0

                    #=
                    if d > d_node[i]
                        d_node[i] = d;
                    end =# 

                    #=
                    #l = ReentrantLock()  
                    lock(l)
                    try
                        #println("lock:", index_th)
                        if d > d_node[i]
                            d_node[i] = d;
                        end        
                    finally
                        #println("unlock:", index_th)
                        unlock(l)
                    end
                    =#

                end # if

            end # for i
            
        end # for j
    
    end  # for c

    #d_node_pre .=  d_node # -> 664 alloc 
    #@inbounds for i = 1 : nNode  # -> 664 alloc 
    #    d_node_pre[i] = d_node[i]
    #end

        #return c_force3, d_node
        #return d_node

    return

end



function cal_contact_force_gpu(c_force3::Array{Float128,2}, CT, instance_pair, cp_index, 
    position::Array{Float64,2}, velo, diag_M, elementMinSize, elementMaxSize, d_max,  
    element_flag, elementmat::Array{Int,2}, bug_report, time_)

    nNode = size(position,2);

    d_lim::Float64 = elementMinSize * 0.3;
    myu = 0.25
    kc_o = 1.0;
    kc_s = 1.0; # self-contact
    Cr_o = 0.0;
    Cr_s = 0.0;

    div3::Float64 = 1.0/3.0


    for c = 1 : length(cp_index)

        cc = cp_index[c]

        #% j -> triangle,   i -> point
        i_instance = instance_pair[c][1]
        j_instance = instance_pair[c][2]
        #println("i_instance:", i_instance, " j_instance:", j_instance)

        c_nodes_i = CT[c].c_nodes_i    # -> 0 alloc
        c_nodes_j = CT[c].c_nodes_j
        c_triangles = CT[c].c_triangles
        c_triangles_eleid = CT[c].c_triangles_eleid
        young = CT[c].young
        nn_i = length(c_nodes_i)
        nn_j = length(c_nodes_j)

        #println("c_nodes_i:", length(c_nodes_i), ", c_triangles:", size(c_triangles,1) )


        #println("c_nodes_i:", c_nodes_i)

        #%--- contact range ---%    
        min_ix = minimum(@view(position[1,c_nodes_i]))
        min_iy = minimum(@view(position[2,c_nodes_i]))
        min_iz = minimum(@view(position[3,c_nodes_i]))
        min_jx = minimum(@view(position[1,c_nodes_j]))
        min_jy = minimum(@view(position[2,c_nodes_j]))
        min_jz = minimum(@view(position[3,c_nodes_j]))

        max_ix = maximum(@view(position[1,c_nodes_i]))
        max_iy = maximum(@view(position[2,c_nodes_i]))
        max_iz = maximum(@view(position[3,c_nodes_i]))
        max_jx = maximum(@view(position[1,c_nodes_j]))
        max_jy = maximum(@view(position[2,c_nodes_j]))
        max_jz = maximum(@view(position[3,c_nodes_j]))

        range_min_x = max(min_ix, min_jx)
        range_min_y = max(min_iy, min_jy)
        range_min_z = max(min_iz, min_jz)
        range_max_x = min(max_ix, max_jx)
        range_max_y = min(max_iy, max_jy)
        range_max_z = min(max_iz, max_jz)

        if range_min_x > range_max_x || range_min_y > range_max_y || range_min_z > range_max_z
            continue
        end

        all_range_min_x = min(min_ix, min_jx)
        all_range_min_y = min(min_iy, min_jy)
        all_range_min_z = min(min_iz, min_jz)
        all_range_max_x = max(max_ix, max_jx)
        all_range_max_y = max(max_iy, max_jy)
        all_range_max_z = max(max_iz, max_jz)
                
        RLx = all_range_max_x - all_range_min_x
        RLy = all_range_max_y - all_range_min_y
        RLz = all_range_max_z - all_range_min_z

        #nn_i = length(c_nodes_i)
        #nn_j = length(c_nodes_j)

        node_map_ix = zeros(Int, nn_i)  
        node_map_iy = zeros(Int, nn_i)  
        node_map_iz = zeros(Int, nn_i)  
        node_map_jx = zeros(Int, nn_j)  
        node_map_jy = zeros(Int, nn_j)  
        node_map_jz = zeros(Int, nn_j)  

        flag_self = false
        kc = kc_o
        Cr = Cr_o
        ddiv = elementMaxSize * 1.1
        if i_instance == j_instance
            flag_self = true
            kc = kc_s
            Cr = Cr_s
            ddiv = elementMaxSize * 0.6
        end
        
        @inbounds for i::Int = 1 : nn_i 
            pos_i_x = position[1, c_nodes_i[i]]
            map_i_x::Int = ceil( (pos_i_x - all_range_min_x) / ddiv )
            node_map_ix[i] = map_i_x

            pos_i_y = position[2, c_nodes_i[i]]
            map_i_y::Int = ceil( (pos_i_y - all_range_min_y) / ddiv )
            node_map_iy[i] = map_i_y

            pos_i_z = position[3, c_nodes_i[i]]
            map_i_z::Int = ceil( (pos_i_z - all_range_min_z) / ddiv )
            node_map_iz[i] = map_i_z
        end
        
        @inbounds for i::Int = 1 : nn_j
            pos_j_x = position[1, c_nodes_j[i]]
            map_j_x::Int = ceil( (pos_j_x - all_range_min_x) / ddiv )
            node_map_jx[i] = map_j_x

            pos_j_y = position[2, c_nodes_j[i]]
            map_j_y::Int = ceil( (pos_j_y - all_range_min_y) / ddiv )
            node_map_jy[i] = map_j_y

            pos_j_z = position[3, c_nodes_j[i]]
            map_j_z::Int = ceil( (pos_j_z - all_range_min_z) / ddiv )
            node_map_jz[i] = map_j_z
        end


        n = size(c_triangles,1)
            #println("n:", n)
        nth = 512 #1024
        nbl = Int32( ceil(n/nth) )
            #println("nth:", nth)
            #println("nbl:", nbl)

        position = cu(position)
        velo = cu(velo)
            
        diag_M = cu(diag_M)
        element_flag = cu(element_flag)
        elementmat = cu(elementmat)
        c_triangles = cu(c_triangles)    
        c_triangles_eleid = cu(c_triangles_eleid)    
        c_nodes_i = cu(c_nodes_i)
        c_nodes_j = cu(c_nodes_j)
        c_force3_gpu = CUDA.zeros(length(velo))
            #println("typeof(position):", typeof(position))
            #println("typeof(velo):", typeof(velo))
            #println("typeof(c_triangles):", typeof(c_triangles))
            #println("typeof(c_nodes_i):", typeof(c_nodes_i))
            #c_force3 = cu(c_force3)
            #println("typeof(c_force3):", typeof(c_force3))
        
        node_map_ix = cu(node_map_ix)
        node_map_iy = cu(node_map_iy)
        node_map_iz = cu(node_map_iz)
        node_map_jx = cu(node_map_jx)
        node_map_jy = cu(node_map_jy)
        node_map_jz = cu(node_map_jz)

        @cuda blocks=nbl threads=nth gpu_contact(n, position, velo, diag_M, element_flag, elementmat, 
                                                c_triangles, c_triangles_eleid, c_nodes_i, c_nodes_j, 
                                                node_map_ix, node_map_iy, node_map_iz, node_map_jx, node_map_jy, node_map_jz, 
                                                range_min_x, range_min_y, range_min_z, range_max_x, range_max_y, range_max_z,
                                                c_force3_gpu, d_lim, young, kc, Cr, myu, flag_self)
        CUDA.synchronize()

        #println("c_force3:", c_force3)

        position = Array(position)
        velo = Array(velo)
        diag_M = Array(diag_M)
        element_flag = Array(element_flag)
        elementmat = Array(elementmat)
        c_triangles = Array(c_triangles)
        c_triangles_eleid = Array(c_triangles_eleid)
        c_nodes_i = Array(c_nodes_i)
        c_nodes_j = Array(c_nodes_j)
        
        c_force3_gpu = Array(c_force3_gpu)
        @inbounds for i = 1 : length(c_force3_gpu)
            c_force3[i] += c_force3_gpu[i]
        end

        
    end  # for c

    return

end


function gpu_contact(n, position, velo, diag_M, element_flag, elementmat, 
                    c_triangles, c_triangles_eleid, c_nodes_i, c_nodes_j,
                    node_map_ix, node_map_iy, node_map_iz, node_map_jx, node_map_jy, node_map_jz, 
                    range_min_x, range_min_y, range_min_z, range_max_x, range_max_y, range_max_z, 
                    c_force3, d_lim, young, kc, Cr, myu, flag_self)

    index::Int64 = (blockIdx().x -1) * blockDim().x + threadIdx().x
    if index > n
        return
    end
    #@cuprint(index, ", ")
    #CUDA.@atomic v[1] += index

    eleid_ = c_triangles_eleid[index]
    if element_flag[ eleid_ ] == 0
        return
    end

    #return

    j0 = c_triangles[index,1]
    j1 = c_triangles[index,2]
    j2 = c_triangles[index,3]

    q0x = position[1,j0]
    q0y = position[2,j0]
    q0z = position[3,j0]
    q1x = position[1,j1]
    q1y = position[2,j1]
    q1z = position[3,j1]
    q2x = position[1,j2]
    q2y = position[2,j2]
    q2z = position[3,j2]

    if q0x < range_min_x && q1x < range_min_x && q2x < range_min_x
        return
    end
    if q0y < range_min_y && q1y < range_min_y && q2y < range_min_y
        return
    end
    if q0z < range_min_z && q1z < range_min_z && q2z < range_min_z
        return
    end

    if q0x > range_max_x && q1x > range_max_x && q2x > range_max_x
        return
    end
    if q0y > range_max_y && q1y > range_max_y && q2y > range_max_y
        return
    end
    if q0z > range_max_z && q1z > range_max_z && q2z > range_max_z
        return
    end

    #=
    cx = (q0x+q1x+q2x) /3.0  #* div3 # /3.0  #
    cy = (q0y+q1y+q2y) /3.0  #* div3 # /3.0
    cz = (q0z+q1z+q2z) /3.0  #* div3 # /3.0
    R0 = gpu3norm(q0x-cx, q0y-cy, q0z-cz)
    R1 = gpu3norm(q1x-cx, q1y-cy, q1z-cz)
    R2 = gpu3norm(q2x-cx, q2y-cy, q2z-cz)
    Rmax = max(max(R0,R1),R2)
    =#

    v1x = q1x - q0x
    v1y = q1y - q0y
    v1z = q1z - q0z

    v2x = q2x - q0x
    v2y = q2y - q0y
    v2z = q2z - q0z

    L1 = gpu3norm(v1x,v1y,v1z);
    L2 = gpu3norm(v2x,v2y,v2z);
    Lmax = max(L1,L2);
    
    nx, ny, nz = gpu3crossNNz(v1x,v1y,v1z, v2x,v2y,v2z)
    #@cuprintln("nx=", nx, ", ny=", ny, ", nz=", nz)
    
    S = 0.0
    d12 = v1x*v2x+v1y*v2y+v1z*v2z
    SS = L1*L1 * L2*L2 - d12*d12
    if SS > 0
        S = 0.5 * sqrt(SS);
    else
        S = 0.5 * L1 * L2
    end
    
    A11 = v1x 
    A21 = v1y 
    A31 = v1z 
    A12 = v2x 
    A22 = v2y 
    A32 = v2z 
    A13 = -nx 
    A23 = -ny 
    A33 = -nz 

    map_j0_x = 1
    map_j0_y = 1
    map_j0_z = 1
    @inbounds for i = 1 : length(c_nodes_j)
        if c_nodes_j[i] == j0
            map_j0_x = node_map_jx[i]
            map_j0_y = node_map_jy[i]
            map_j0_z = node_map_jz[i]
            break
        end
    end 

    e1_ = elementmat[1, eleid_]
    e2_ = elementmat[2, eleid_]
    e3_ = elementmat[3, eleid_]
    e4_ = elementmat[4, eleid_]
    e5_ = elementmat[5, eleid_]
    e6_ = elementmat[6, eleid_]
    e7_ = elementmat[7, eleid_]
    e8_ = elementmat[8, eleid_]


    @inbounds for  k = 1 : length(c_nodes_i)

        i = c_nodes_i[k]

        if (abs(map_j0_x - node_map_ix[k]) > 1  ||
            abs(map_j0_y - node_map_iy[k]) > 1  ||
            abs(map_j0_z - node_map_iz[k]) > 1  )
            
            continue
        end 

        if flag_self
            if ( i == e1_ ||
                i == e2_ ||
                i == e3_ ||
                i == e4_ ||
                i == e5_ ||
                i == e6_ ||
                i == e7_ ||
                i == e8_  )
                continue
            end 
        end

        px = position[1,i]
        py = position[2,i]
        pz = position[3,i]

        #=
        if px < range_min_x || py < range_min_y || pz < range_min_z   # -> slow?
            continue
        end
        if px > range_max_x || py > range_max_y || pz > range_max_z
            continue
        end
        =#

        #dpc = gpu3norm(px-cx,py-cy,pz-cz)  # -> slow?
        #if dpc >= Rmax
            #continue  
        #end   

        bx = px - q0x
        by = py - q0y
        bz = pz - q0z

        x1, x2, d = gpu3SolveAb(A11,A21,A31, A12,A22,A32, A13,A23,A33,
                                    bx,by,bz)

        if (0.0 <= x1 && 0.0 <= x2 && x1 + x2 <= 1.0 &&
            d > 0.0 && d <= d_lim)

            #if d - d_node_pre[i] > d_max
            #    d = d_node_pre[i] + d_max;
            #end

            k = young * S / Lmax  * kc
            F = k * d

            fx = F * nx
            fy = F * ny
            fz = F * nz
            #@cuprintln(index, ", fx=", fx, ", fy=", fy, ", fz=", fz)
                    
            #@cuprintln("typeof(vex)=", typeof(vex))
            #@cuprintln("typeof(vey)=", typeof(vey))
            #@cuprintln("typeof(vez)=", typeof(vez))
            #@cuprintln("typeof(nx)=", typeof(nx))
            #@cuprintln("typeof(ny)=", typeof(ny))
            #@cuprintln("typeof(nz)=", typeof(nz))

            #@cuprintln(index, ", d=", d)
            vx = velo[i*3-2] - velo[j0*3-2]
            vy = velo[i*3-1] - velo[j0*3-1]
            vz = velo[i*3] - velo[j0*3]
            mag_v = gpu3norm(vx,vy,vz)
            
            vex = 0.0
            vey = 0.0
            vez = 0.0
            dot_ve_n = 0.0
                    
            if mag_v > 0.0
                vex = vx / mag_v
                vey = vy / mag_v
                vez = vz / mag_v  
                #@cuprintln("vex=", vex, ", vey=", vey, ", vez=", vez)
                dot_ve_n = vex*nx + vey*ny + vez*nz  
            end
            #@cuprintln("vex=", vex)
            #@cuprintln("vey=", vey)
            #@cuprintln("vez=", vez)
            #@cuprintln("vex=", vex, ", vey=", vey, ", vez=", vez)

            C = 2 * sqrt( diag_M[i] * k ) * Cr;
            fc_x = -C * vx
            fc_y = -C * vy
            fc_z = -C * vz
            
            vsx = vex - dot_ve_n * nx
            vsy = vey - dot_ve_n * ny
            vsz = vez - dot_ve_n * nz
            
            fric_x = -myu * F * vsx
            fric_y = -myu * F * vsy
            fric_z = -myu * F * vsz
            fx += fric_x + fc_x
            fy += fric_y + fc_y
            fz += fric_z + fc_z

            CUDA.@atomic  c_force3[1+(i-1)*3] +=  fx
            CUDA.@atomic  c_force3[2+(i-1)*3] +=  fy
            CUDA.@atomic  c_force3[3+(i-1)*3] +=  fz

            CUDA.@atomic  c_force3[1+(j0-1)*3] +=  -fx  / 3.0
            CUDA.@atomic  c_force3[2+(j0-1)*3] +=  -fy  / 3.0
            CUDA.@atomic  c_force3[3+(j0-1)*3] +=  -fz  / 3.0

            CUDA.@atomic  c_force3[1+(j1-1)*3] +=  -fx  / 3.0
            CUDA.@atomic  c_force3[2+(j1-1)*3] +=  -fy  / 3.0
            CUDA.@atomic  c_force3[3+(j1-1)*3] +=  -fz  / 3.0

            CUDA.@atomic  c_force3[1+(j2-1)*3] +=  -fx  / 3.0
            CUDA.@atomic  c_force3[2+(j2-1)*3] +=  -fy  / 3.0
            CUDA.@atomic  c_force3[3+(j2-1)*3] +=  -fz  / 3.0

            #if d > d_node[i]
                #CUDA.@atomic  d_node[i] += d;
                #CUDA.atomic_add!(pointer(d_node[i]), d)
            #    d_node[i] = d;  # does not work
            #end

        end # if

    end#for k 

    return     

end




function my3norm(v1)   # @inline -> even
    v = sqrt(v1[1]*v1[1]+v1[2]*v1[2]+v1[3]*v1[3])
    return v
end

function my3norm(b1::Float64,b2::Float64,b3::Float64)   # @inline -> even
    v::Float64 = sqrt(b1*b1+b2*b2+b3*b3)
    return v
end

function gpu3norm(b1,b2,b3)
    v = sqrt(b1*b1+b2*b2+b3*b3)
    return v
end


function my3cross(a,b)
    n = zeros(3);  #@MVector 
    n[1] = a[2]*b[3] - a[3]*b[2];
    n[2] = a[3]*b[1] - a[1]*b[3];
    n[3] = a[1]*b[2] - a[2]*b[1];
    return n
end

function my3crossN(a,b,n)   # reduce alloc and memory compared to my3cross(a,b)
    #n = zeros(3);  #@MVector 
    n[1] = a[2]*b[3] - a[3]*b[2];
    n[2] = a[3]*b[1] - a[1]*b[3];
    n[3] = a[1]*b[2] - a[2]*b[1];
    #return n
    return
end

function my3crossNNz(a,b,n)   # Normalize
    
    n1 = a[2]*b[3] - a[3]*b[2];
    n2 = a[3]*b[1] - a[1]*b[3];
    n3 = a[1]*b[2] - a[2]*b[1];
    mag_n = sqrt( n1*n1 + n2*n2 + n3*n3 )
    n[1] = n1 / mag_n
    n[2] = n2 / mag_n
    n[3] = n3 / mag_n
    
    return
end

#function my3crossNNz(a1,a2,a3, b1,b2,b3)   # Normalize
function my3crossNNz(a1::Float64, a2::Float64, a3::Float64,  b1::Float64, b2::Float64, b3::Float64)    
    
    n1 = a2*b3 - a3*b2
    n2 = a3*b1 - a1*b3
    n3 = a1*b2 - a2*b1
    mag_n = sqrt( n1*n1 + n2*n2 + n3*n3 )
    n1 = n1 / mag_n
    n2 = n2 / mag_n
    n3 = n3 / mag_n
    
    return n1, n2, n3
end

function gpu3crossNNz(a1, a2, a3,  b1, b2, b3)    
    
    n1 = a2*b3 - a3*b2
    n2 = a3*b1 - a1*b3
    n3 = a1*b2 - a2*b1
    mag_n = sqrt( n1*n1 + n2*n2 + n3*n3 )
    n1 = n1 / mag_n
    n2 = n2 / mag_n
    n3 = n3 / mag_n
    
    return n1, n2, n3
end

function my3det(m3)
    v = ( m3[1,1]*m3[2,2]*m3[3,3] 
        + m3[1,2]*m3[2,3]*m3[3,1] 
        + m3[1,3]*m3[2,1]*m3[3,2] 
        - m3[1,1]*m3[2,3]*m3[3,2] 
        - m3[1,2]*m3[2,1]*m3[3,3] 
        - m3[1,3]*m3[2,2]*m3[3,1] )
    return v
end

function my3inv(m3)
    v = ( m3[1,1]*m3[2,2]*m3[3,3] 
        + m3[1,2]*m3[2,3]*m3[3,1] 
        + m3[1,3]*m3[2,1]*m3[3,2] 
        - m3[1,1]*m3[2,3]*m3[3,2] 
        - m3[1,2]*m3[2,1]*m3[3,3] 
        - m3[1,3]*m3[2,2]*m3[3,1] )

    im = zeros(3,3); #@MMatrix 
    #im = Matrix{Float64}(undef, 3, 3) # -> same
    im[1,1] = ( m3[2,2]*m3[3,3] - m3[2,3]*m3[3,2] ) / v;
    im[2,1] = ( m3[2,3]*m3[3,1] - m3[2,1]*m3[3,3] ) / v;
    im[3,1] = ( m3[2,1]*m3[3,2] - m3[2,2]*m3[3,1] ) / v;

    im[1,2] = ( m3[1,3]*m3[3,2] - m3[1,2]*m3[3,3] ) / v;
    im[2,2] = ( m3[1,1]*m3[3,3] - m3[1,3]*m3[3,1] ) / v;
    im[3,2] = ( m3[1,2]*m3[3,1] - m3[1,1]*m3[3,2] ) / v;

    im[1,3] = ( m3[1,2]*m3[2,3] - m3[1,3]*m3[2,2] ) / v;
    im[2,3] = ( m3[1,3]*m3[2,1] - m3[1,1]*m3[2,3] ) / v;
    im[3,3] = ( m3[1,1]*m3[2,2] - m3[1,2]*m3[2,1] ) / v;

    return im;
end

function my3SolveAb(A,b,x)
    v = ( A[1,1]*A[2,2]*A[3,3] 
        + A[1,2]*A[2,3]*A[3,1] 
        + A[1,3]*A[2,1]*A[3,2] 
        - A[1,1]*A[2,3]*A[3,2] 
        - A[1,2]*A[2,1]*A[3,3] 
        - A[1,3]*A[2,2]*A[3,1] )

    #im = Matrix{Float64}(undef, 3, 3) # -> same
#=    im = zeros(3,3); #@MMatrix 
    im[1,1] = ( A[2,2]*A[3,3] - A[2,3]*A[3,2] ) / v;
    im[2,1] = ( A[2,3]*A[3,1] - A[2,1]*A[3,3] ) / v;
    im[3,1] = ( A[2,1]*A[3,2] - A[2,2]*A[3,1] ) / v;

    im[1,2] = ( A[1,3]*A[3,2] - A[1,2]*A[3,3] ) / v;
    im[2,2] = ( A[1,1]*A[3,3] - A[1,3]*A[3,1] ) / v;
    im[3,2] = ( A[1,2]*A[3,1] - A[1,1]*A[3,2] ) / v;

    im[1,3] = ( A[1,2]*A[2,3] - A[1,3]*A[2,2] ) / v;
    im[2,3] = ( A[1,3]*A[2,1] - A[1,1]*A[2,3] ) / v;
    im[3,3] = ( A[1,1]*A[2,2] - A[1,2]*A[2,1] ) / v;

    x[1] = im[1,1]*b[1] + im[1,2]*b[2] + im[1,3]*b[3]
    x[2] = im[2,1]*b[1] + im[2,2]*b[2] + im[2,3]*b[3]
    x[3] = im[3,1]*b[1] + im[3,2]*b[2] + im[3,3]*b[3]
=#
    im11 = ( A[2,2]*A[3,3] - A[2,3]*A[3,2] ) / v;
    im21 = ( A[2,3]*A[3,1] - A[2,1]*A[3,3] ) / v;
    im31 = ( A[2,1]*A[3,2] - A[2,2]*A[3,1] ) / v;

    im12 = ( A[1,3]*A[3,2] - A[1,2]*A[3,3] ) / v;
    im22 = ( A[1,1]*A[3,3] - A[1,3]*A[3,1] ) / v;
    im32 = ( A[1,2]*A[3,1] - A[1,1]*A[3,2] ) / v;

    im13 = ( A[1,2]*A[2,3] - A[1,3]*A[2,2] ) / v;
    im23 = ( A[1,3]*A[2,1] - A[1,1]*A[2,3] ) / v;
    im33 = ( A[1,1]*A[2,2] - A[1,2]*A[2,1] ) / v;

    x[1] = im11*b[1] + im12*b[2] + im13*b[3]
    x[2] = im21*b[1] + im22*b[2] + im23*b[3]
    x[3] = im31*b[1] + im32*b[2] + im33*b[3]

    return 
end

function my3SolveAb(A, bx::Float64, by::Float64, bz::Float64, x)
    v = ( A[1,1]*A[2,2]*A[3,3] 
        + A[1,2]*A[2,3]*A[3,1] 
        + A[1,3]*A[2,1]*A[3,2] 
        - A[1,1]*A[2,3]*A[3,2] 
        - A[1,2]*A[2,1]*A[3,3] 
        - A[1,3]*A[2,2]*A[3,1] )

    im11 = ( A[2,2]*A[3,3] - A[2,3]*A[3,2] ) / v;
    im21 = ( A[2,3]*A[3,1] - A[2,1]*A[3,3] ) / v;
    im31 = ( A[2,1]*A[3,2] - A[2,2]*A[3,1] ) / v;

    im12 = ( A[1,3]*A[3,2] - A[1,2]*A[3,3] ) / v;
    im22 = ( A[1,1]*A[3,3] - A[1,3]*A[3,1] ) / v;
    im32 = ( A[1,2]*A[3,1] - A[1,1]*A[3,2] ) / v;

    im13 = ( A[1,2]*A[2,3] - A[1,3]*A[2,2] ) / v;
    im23 = ( A[1,3]*A[2,1] - A[1,1]*A[2,3] ) / v;
    im33 = ( A[1,1]*A[2,2] - A[1,2]*A[2,1] ) / v;

    x[1] = im11*bx + im12*by + im13*bz
    x[2] = im21*bx + im22*by + im23*bz
    x[3] = im31*bx + im32*by + im33*bz

    return 
end

function my3SolveAb(A11::Float64, A21::Float64, A31::Float64, 
    A12::Float64, A22::Float64, A32::Float64, 
    A13::Float64, A23::Float64, A33::Float64, 
    bx::Float64, by::Float64, bz::Float64)

    v = ( A11*A22*A33 
        + A12*A23*A31 
        + A13*A21*A32 
        - A11*A23*A32 
        - A12*A21*A33 
        - A13*A22*A31 )

    #im11 = ( A22*A33 - A23*A32 ) / v;
    im11 = A22*A33 - A23*A32 
    im21 = A23*A31 - A21*A33 
    im31 = A21*A32 - A22*A31 

    im12 = A13*A32 - A12*A33
    im22 = A11*A33 - A13*A31
    im32 = A12*A31 - A11*A32

    im13 = A12*A23 - A13*A22
    im23 = A13*A21 - A11*A23
    im33 = A11*A22 - A12*A21

    #x1 = im11*bx + im12*by + im13*bz
    x1 = (im11*bx + im12*by + im13*bz) / v
    x2 = (im21*bx + im22*by + im23*bz) / v
    x3 = (im31*bx + im32*by + im33*bz) / v

    return x1, x2, x3 
end

function gpu3SolveAb(A11, A21, A31, 
                    A12, A22, A32, 
                    A13, A23, A33, 
                    bx, by, bz)

    v = ( A11*A22*A33 
        + A12*A23*A31 
        + A13*A21*A32 
        - A11*A23*A32 
        - A12*A21*A33 
        - A13*A22*A31 )

    im11 = A22*A33 - A23*A32 
    im21 = A23*A31 - A21*A33 
    im31 = A21*A32 - A22*A31 

    im12 = A13*A32 - A12*A33
    im22 = A11*A33 - A13*A31
    im32 = A12*A31 - A11*A32

    im13 = A12*A23 - A13*A22
    im23 = A13*A21 - A11*A23
    im33 = A11*A22 - A12*A21

    x1 = (im11*bx + im12*by + im13*bz) / v
    x2 = (im21*bx + im22*by + im23*bz) / v
    x3 = (im31*bx + im32*by + im33*bz) / v

    return x1, x2, x3 
end


#function cal_node_stress_strain(nNode, elementmat, integ_num, integ_stress, integ_strain)
function cal_node_stress_strain(nNode, elementmat, integ_num, integ_data)

    integ_stress = integ_data.integ_stress'; # row major
    integ_strain = integ_data.integ_strain';
    integ_plastic_strain = integ_data.integ_plastic_strain;
    integ_eq_plastic_strain = integ_data.integ_eq_plastic_strain;
    integ_triax_stress = integ_data.integ_triax_stress;

    node_stress = zeros(nNode, 6);
    node_strain = zeros(nNode, 6);
    node_plastic_strain = zeros(nNode, 6);
    node_eq_plastic_strain = zeros(nNode);
    node_triax_stress = zeros(nNode);

    nElement = size(elementmat,2); #size(elementmat,1)
    element_stress = zeros(nElement, 6);
    element_strain = zeros(nElement, 6);
    element_eq_plastic_strain = zeros(nElement);
    element_triax_stress = zeros(nElement);

    for i = 1 : nElement
        if integ_num == 1
            element_stress[i,:] = integ_stress[i,:];
            element_strain[i,:] = integ_strain[i,:];
            element_eq_plastic_strain[i] = integ_eq_plastic_strain[i];
            element_triax_stress[i] = integ_triax_stress[i];
        else
            element_stress[i,:] = sum!(zeros(1,6), integ_stress[(1:integ_num).+(i-1)*integ_num, :]  ) / integ_num;
            element_strain[i,:] = sum!(zeros(1,6), integ_strain[(1:integ_num).+(i-1)*integ_num, :]  ) / integ_num;
            element_eq_plastic_strain[i] = sum( integ_eq_plastic_strain[(1:integ_num).+(i-1)*integ_num]  ) / integ_num;
            element_triax_stress[i] = sum( integ_triax_stress[(1:integ_num).+(i-1)*integ_num]  ) / integ_num;
        end
    end

    for i = 1 : nElement
        for k = 1 : 8
            #node_stress[elementmat[i,k],:] += element_stress[i, :]
            #node_strain[elementmat[i,k],:] += element_strain[i, :]
            #node_eq_plastic_strain[elementmat[i,k]] += element_eq_plastic_strain[i]
            #node_triax_stress[elementmat[i,k]] += element_triax_stress[i]
            node_stress[elementmat[k,i],:] += element_stress[i, :]
            node_strain[elementmat[k,i],:] += element_strain[i, :]
            node_eq_plastic_strain[elementmat[k,i]] += element_eq_plastic_strain[i]
            node_triax_stress[elementmat[k,i]] += element_triax_stress[i]
        end
    end

    inc_num = zeros(nNode, 1);
    for i = 1 : nElement
       #inc_num[elementmat[i,:]] .+= 1 
       inc_num[elementmat[:,i]] .+= 1 
    end
    #println("inc_num:", inc_num)

    for i = 1 : nNode
        node_stress[i,:] ./= inc_num[i]
        node_strain[i,:] ./= inc_num[i]
        node_eq_plastic_strain[i] /= inc_num[i]
        node_triax_stress[i] /= inc_num[i]
    end

    node_mises_stress = zeros(nNode);
    for i = 1 : nNode
        ox = node_stress[i,1];
        oy = node_stress[i,2];
        oz = node_stress[i,3];
        txy = node_stress[i,4];
        tyz = node_stress[i,5];
        txz = node_stress[i,6];
        node_mises_stress[i] = sqrt( 0.5*( (ox-oy)^2 + (oy-oz)^2 + (ox-oz)^2 + 6 * (txy^2 + tyz^2 + txz^2) ) );
    end

    #return node_stress, node_strain, node_mises_stress

    node_data = NodeDataType(node_stress, node_strain, node_plastic_strain, node_eq_plastic_strain,
                             node_mises_stress, node_triax_stress)

    return node_data
end

#=
function drawElement(coordmat, elementmat)

    xp = [0, 0, 0, 0, 1, 1, 1, 1]
    yp = [0, 1, 0, 1, 0, 0, 1, 1]
    zp = [0, 0, 1, 1, 1, 0, 0, 1]
    connections = [(1,2,3), (4,2,3), (4,7,8), (7,5,6), (2,4,7), (1,6,2), (2,7,6), (7,8,5), (4,8,5), (4,5,3), (1,6,3), (6,3,5)]
    
    xe = [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0]
    ye = [0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1]
    ze = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1]
    
    plot(xe,ye,ze; lc=:black, lw=0.5, lims=(-0.25,1.25))
    scatter!((0.5,0.5,0.5); c=:red, ms=6, msw=0.1)
    mesh3d!(xp,yp,zp; connections, proj_type=:persp, fc=:lime, lc=:lime, fa=0.3, lw=0)
    gui()

end =#


function drawGraph(output)

    plot(output)
    gui()

end


#function write_vtk(index, coordmat, elementmat, element_flag, disp, node_data)
function write_vtk(index, coordmat, elementmat, element_flag, disp, velo, node_data)

    node_stress = node_data.node_stress;
    node_strain = node_data.node_strain;
    node_mises_stress = node_data.node_mises_stress;
    node_eq_plastic_strain = node_data.node_eq_plastic_strain;
    node_triax_stress = node_data.node_triax_stress;

    nNode = size(coordmat,2)  #size(coordmat,1)
    nElement = size(elementmat,2) #size(elementmat,1)
    disp3 = reshape(disp,3,nNode)'
    velo3 = reshape(velo,3,nNode)'

    for i = 1 : nNode
        for j = 1 : 3
            if abs(disp3[i,j]) < 1E-16
                disp3[i,j] = 0;
            end
            if abs(velo3[i,j]) < 1E-16
                velo3[i,j] = 0;
            end
        end

        for j = 1 : 6
            if abs(node_stress[i,j]) < 1E-16
                node_stress[i,j] = 0;
            end
            if abs(node_strain[i,j]) < 1E-16
                node_strain[i,j] = 0;
            end
        end

        if abs(node_mises_stress[i]) < 1E-16
            node_mises_stress[i] = 0;
        end
        if abs(node_eq_plastic_strain[i]) < 1E-16
            node_eq_plastic_strain[i] = 0;
        end
        if abs(node_triax_stress[i]) < 1E-16
            node_triax_stress[i] = 0;
        end
    end

    if !ispath("temp")
        mkdir("temp")
    end

    fname = @sprintf("temp\\file%03d.vtk", index)
    out = open(fname,"w")
    println(out,"# vtk DataFile Version 2.0")
    println(out,"Test")
    println(out,"ASCII")
    println(out,"DATASET UNSTRUCTURED_GRID")
    println(out,"POINTS ", nNode, " float")
    
    for i = 1 : nNode
        #println(out, @sprintf("%1.6e %1.6e %1.6e", coordmat[i,1], coordmat[i,2], coordmat[i,3]) )
        println(out, @sprintf("%1.6e %1.6e %1.6e", coordmat[1,i], coordmat[2,i], coordmat[3,i]) )
    end

    draw_element_num = sum(element_flag)

    #println(out,"CELLS ", nElement, " ", nElement*(8+1))
    println(out,"CELLS ", draw_element_num, " ", draw_element_num*(8+1))
    e = elementmat;
    e = e .- 1;
    for i = 1 : nElement
        if element_flag[i] == 1
            println(out, @sprintf("8 %d %d %d %d %d %d %d %d", 
            e[1,i], e[2,i], e[3,i], e[4,i], e[5,i], e[6,i], e[7,i], e[8,i]) )
            #e[i,1], e[i,2], e[i,3], e[i,4], e[i,5], e[i,6], e[i,7], e[i,8]) )
        end
    end

    #println(out,"CELL_TYPES ", nElement)
    println(out,"CELL_TYPES ", draw_element_num)
    for i = 1 : draw_element_num  #nElement
        println(out, 12)
    end

    
    println(out,"POINT_DATA ", nNode)
    println(out,"VECTORS DISPLACEMENT float")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e %1.6e %1.6e", disp3[i,1], disp3[i,2], disp3[i,3]) )
    end


    println(out,"SCALARS Vx float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", velo3[i,1] ) )
    end

    println(out,"SCALARS Vy float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", velo3[i,2] ) )
    end

    println(out,"SCALARS Vz float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", velo3[i,3] ) )
    end


    println(out,"SCALARS E11 float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_strain[i,1] ) )
    end

    println(out,"SCALARS E22 float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_strain[i,2] ) )
    end

    println(out,"SCALARS E33 float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_strain[i,3] ) )
    end

    println(out,"SCALARS E12 float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_strain[i,4] ) )
    end

    println(out,"SCALARS E23 float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_strain[i,5] ) )
    end

    println(out,"SCALARS E13 float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_strain[i,6] ) )
    end

    println(out,"SCALARS EQ_PSTRAIN float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_eq_plastic_strain[i] ) )
    end


    println(out,"SCALARS S11 float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_stress[i,1] ) )
    end

    println(out,"SCALARS S22 float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_stress[i,2] ) )
    end

    println(out,"SCALARS S33 float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_stress[i,3] ) )
    end

    println(out,"SCALARS S12 float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_stress[i,4] ) )
    end

    println(out,"SCALARS S23 float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_stress[i,5] ) )
    end

    println(out,"SCALARS S13 float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_stress[i,6] ) )
    end

    println(out,"SCALARS MISES_STRESS float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_mises_stress[i] ) )
    end

    println(out,"SCALARS TRIAX_STRESS float 1")
    println(out,"LOOKUP_TABLE default")
    for i = 1 : nNode
        println(out, @sprintf("%1.6e", node_triax_stress[i] ) )
    end

    close(out)

end


#gr()


#println(ARGS)
#println(ARGS[1])
#@time res = hakai(ARGS[1])

#@time res = hakai("C:\\Users\\under\\matlab\\HAKAI\\input\\Tensile5e.inp")

function main()
    #@time res = hakai("C:\\Users\\under\\matlab\\HAKAI\\input\\Tensile5e.inp")
    #@time 
    res = hakai(ARGS[1])
end

main()
