function [ A,b,u,freenode] = fractures_matrix( node,fractures,intersections,materials,fracure_lengths)
%FRACTURES_MATRIX Summary of this function goes here
%   Detailed explanation goes here
tmp_len=[0 cumsum(fracure_lengths)];
n=length(fractures);
tmp_matrices=cell(n,3);

for i=1:n
    [ tmp_matrices{i,1},tmp_matrices{i,2}] = FEM1D_mass(node,fractures{i});
    [ tmp_matrices{i,3},b,u,freenode ] = FEM1D_2(node,fractures{i});
end

[ G ] = intersects_glue_part( intersections,materials,fracure_lengths );

A=G;

for i=1:n
    A(tmp_len(i)+1:tmp_len(i+1),tmp_len(i)+1:tmp_len(i+1))=A(tmp_len(i)+1:tmp_len(i+1),tmp_len(i)+1:tmp_len(i+1))+tmp_matrices{i,1}+tmp_matrices{i,2}+tmp_matrices{i,3};
end

b=zeros(length(A),1);
u=zeros(length(A),1);
freenode=(1:length(A))';
end

