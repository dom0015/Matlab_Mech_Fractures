function [ F_stif,b,u,freenode,F_mass] = fractures_matrix_modif( node,fractures,intersections,materials,fracure_lengths)
%FRACTURES_MATRIX Summary of this function goes here
%   Detailed explanation goes here
tmp_len=[0 cumsum(fracure_lengths)];
n=length(fractures);
tmp_matrices=cell(n,3);

for i=1:n
    [ tmp_matrices{i,1},tmp_matrices{i,2}] = a_hyd.FEM1D_mass(node,fractures{i});
    [ tmp_matrices{i,3},b,u,freenode ] = a_hyd.FEM1D_2(node,fractures{i});
end

F_mass = a_hyd.intersects_glue_part( intersections,materials,fracure_lengths );
F_stif = F_mass;

for i=1:n
    F_mass(tmp_len(i)+1:tmp_len(i+1),tmp_len(i)+1:tmp_len(i+1))=...
        F_mass(tmp_len(i)+1:tmp_len(i+1),tmp_len(i)+1:tmp_len(i+1))+tmp_matrices{i,1}+tmp_matrices{i,2};
    F_stif(tmp_len(i)+1:tmp_len(i+1),tmp_len(i)+1:tmp_len(i+1))=...
        F_stif(tmp_len(i)+1:tmp_len(i+1),tmp_len(i)+1:tmp_len(i+1))+tmp_matrices{i,3};
end

b=zeros(length(F_mass),1);
u=zeros(length(F_mass),1);
freenode=(1:length(F_mass))';
end

