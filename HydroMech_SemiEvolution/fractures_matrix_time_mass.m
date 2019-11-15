function [freenode,F_mass] = fractures_matrix_time_mass( node,fractures,intersections,materials,fracure_lengths)
%FRACTURES_MATRIX Summary of this function goes here
%   Detailed explanation goes here
tmp_len=[0 cumsum(fracure_lengths)];
n=length(fractures);
tmp_matrices=cell(n,3);

for i=1:n
     [ tmp_matrices{i,1}] = a_hyd.FEM1D_time(node,fractures{i});
end

F_mass = a_hyd.intersects_glue_part( intersections,materials,fracure_lengths )*0;

for i=1:n
    F_mass(tmp_len(i)+1:tmp_len(i+1),tmp_len(i)+1:tmp_len(i+1))=...
        F_mass(tmp_len(i)+1:tmp_len(i+1),tmp_len(i)+1:tmp_len(i+1))+tmp_matrices{i,1};
end
freenode=(1:length(F_mass))';
end

