function [x] = solve_penalized_operator(A,B,rho,b)
%SOLVE_PENALIZED_OPERATOR Summary of this function goes here
%   Detailed explanation goes here
[n,m]=size(B);
r_s=sqrt(rho);
M=[A B';B -speye(n,n)/rho];
f=[b;zeros(n,1)];
x=M\f;
x=x(1:m);
end

