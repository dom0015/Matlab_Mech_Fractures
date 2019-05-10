function [ B, freeNode, b, u ] = matrices_assembling( A, I, M, F, freeNode, freeNode_f, b, b_f, u0, u_f )
%MATRICES_ASSEMBLING Summary of this function goes here
%   Detailed explanation goes here
N = length(A);
freeNode = [freeNode; freeNode_f+N ];
b = [b; b_f];
u = [u0; u_f];
N_i = length(F) - size(I,1);
Z = zeros(N_i,N);
I = [I; Z];
B = [A+M -I'; -I F];

end

