
% rel=1.0e-14;
% rho0=1;
% betarho=1.001;
% Gama = 1;
% maxiter_cg = 10000;
% M_start=0.1;
% type='M';


rhoratio=eigs(G'*G,1)/eigs(F,1);

rel=1.0e-10;
rho0=0.01;
betarho=2;
Gama = 1;
maxiter_cg = 3000;
M_start=100;
type='m';


tic
[lambda_ker] = SMALSE(F,d-F*lambda_ImGt,G,c_ker,idx_no_bounds,rel,rho0,betarho,Gama,M_start,maxiter_cg,type);
toc