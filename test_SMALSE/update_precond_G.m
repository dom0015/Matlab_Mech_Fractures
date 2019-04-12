function [G,Gt,Q] = update_precond_G(u,recalculate_B,C_inv,Ct_inv,n)
%UPDATE_PRECOND_G Summary of this function goes here
%   Detailed explanation goes here
global u_p
%B=recalculate_B(C_inv(u(1:n)));
B=recalculate_B(u_p);
[B,R] = my_QR(B);
R=pinv(full(R));
BtB=B'*B;
G=@(x)B*C_inv(x(1:n))+R*x((n+1):end);
Q=@(x)[Ct_inv(BtB*C_inv(x(1:n)))+Ct_inv(B'*R*x((n+1):end));R'*B*C_inv(x(1:n))+R'*R*x((n+1):end)];
Gt=@(x)[Ct_inv(B'*x);R'*x];
end