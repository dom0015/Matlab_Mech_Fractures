function [problem_setting] = assembly_FETI_frac_rhs(problem_setting,frac_press)
%ASSEMBLY_FETI_FRAC_RHS Summary of this function goes here
%   Detailed explanation goes here
sub_nodes=problem_setting.sub_nodes;
sub_elem=problem_setting.sub_elem;
fracture_elem_map=problem_setting.fracture_elem_map;
divide_neumann_boundary=problem_setting.divide_neumann_boundary;
Neumann_normalx=problem_setting.Neumann_normalx;
Neumann_normaly=problem_setting.Neumann_normaly;
Neumann_boundaries=problem_setting.N_bound_fractures;

sumbdomains_FETI=length(sub_nodes);
N_bound_value=cell(2,1);
tmp_nx=Neumann_normalx;
tmp_ny=Neumann_normaly;

for i=1:length(fracture_elem_map)
    [tmp_nx,tmp_ny] = fracture_boundary_values(fracture_elem_map{i},tmp_nx,tmp_ny,frac_press{i},frac_press{i});
end
N_bound_value{1}=tmp_nx;
N_bound_value{2}=tmp_ny;

[sub_boundary,sub_boundary_val]=divide_neumann_boundary(Neumann_boundaries,N_bound_value);

%% Elasticity sub-Matrices assembly ---------------------------------------
b=cell(sumbdomains_FETI,1);
for i=1:sumbdomains_FETI
    [b{i}]=elasticity_assembly_neumann_only(sub_nodes{i},sub_elem{i},sub_boundary{i},sub_boundary_val{i});
end
b=cat(1,b{:});
problem_setting.b_frac=b;
end