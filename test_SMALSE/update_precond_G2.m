function [G,Gt,Q] = update_precond_G2(u,recalculate_B,C_inv,Ct_inv,n)
%UPDATE_PRECOND_G Summary of this function goes here
%   Detailed explanation goes here
global u_p
tmp=C_inv(u);
B=recalculate_B(u_p);
GG=[B eye(length(u)-n)];
QQ=GG'*GG;
G=@(x)GG*C_inv(x);
Q=@(x)Ct_inv(QQ*C_inv(x));
Gt=@(x)Ct_inv(GG'*x);
end

