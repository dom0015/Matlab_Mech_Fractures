function [k] = proj_B_larger(u,B,c,rho0)
%PROJ_B_LARGER Summary of this function goes here
%   Detailed explanation goes here
rho0=rho0*ones(size(B,1),1);
idx=B*u-c<rho0;
A=B(idx,:)'*B(idx,:);
k=A\(-A*u+B(idx,:)'*(c(idx)+rho0(idx)));


x = quadprog(speye(size(B,2)),0*u,-B,B*u-c-rho0);

end

