function [L] = chol_reg(A,tol)
%CHOL_REG Summary of this function goes here
%   Detailed explanation goes here
A=(A+A')/2;
I=speye(size(A,1));
A_reg=A+tol*I;
L=chol(A_reg);
end

