function [problem_setting] = assembly_FETI_frac_rhs(problem_setting,frac_press,minus_ugrad)
%ASSEMBLY_FETI_FRAC_RHS Summary of this function goes here
%   Detailed explanation goes here
% figure(204); plot(cell2mat(frac_press)); hold on
% figure(205); plot(minus_ugrad); hold on
sub_nodes=problem_setting.sub_nodes;
sub_elem=problem_setting.sub_elem;
fracture_elem_map=problem_setting.fracture_elem_map;
divide_neumann_boundary=problem_setting.divide_neumann_boundary;
divide_neumann_boundary_and_force=problem_setting.divide_neumann_boundary_and_force;
Neumann_normalx=problem_setting.Neumann_normalx;
Neumann_normaly=problem_setting.Neumann_normaly;
Neumann_boundaries=problem_setting.N_bound_fractures;

sumbdomains_FETI=length(sub_nodes);
N_bound_value=cell(2,1);
tmp_nx=Neumann_normalx;
tmp_ny=Neumann_normaly;

if size(frac_press,2)==1
    for i=1:length(fracture_elem_map)
        [tmp_nx,tmp_ny] = fracture_boundary_values(fracture_elem_map{i},tmp_nx,tmp_ny,frac_press{i},frac_press{i});
    end
else
    for i=1:length(fracture_elem_map)
        [tmp_nx,tmp_ny] = feti.fracture_boundary_values(fracture_elem_map{i},tmp_nx,tmp_ny,frac_press{i,1},frac_press{i,2});
    end
end
N_bound_value{1}=tmp_nx;
N_bound_value{2}=tmp_ny;

if isempty(minus_ugrad)
    [sub_boundary,sub_boundary_val]=divide_neumann_boundary(Neumann_boundaries,N_bound_value);
    %% Elasticity sub-Matrices assembly ---------------------------------------
    b=cell(sumbdomains_FETI,1);
    for i=1:sumbdomains_FETI
        [b{i}]=elasticity_assembly_neumann_only(sub_nodes{i},sub_elem{i},sub_boundary{i},sub_boundary_val{i});
    end
else
    [sub_boundary,sub_boundary_val,volume_forces]=divide_neumann_boundary_and_force(Neumann_boundaries,N_bound_value,minus_ugrad);
    %% Elasticity sub-Matrices assembly ---------------------------------------
    b=cell(sumbdomains_FETI,1);
    for i=1:sumbdomains_FETI
        [b{i}]=a_ela.elasticity_assembly_neumann_and_force_only(sub_nodes{i},sub_elem{i},sub_boundary{i},sub_boundary_val{i},volume_forces{i});
    end
    
end
b=cat(1,b{:});
problem_setting.b_frac=b;
end