function [D,elast_problem,x_elast,ncg] = SMALSE_solver_new(elast_problem,SMALSE_params)
%SMALSE_SOLVER Summary of this function goes here
%   Detailed explanation goes here
b=(elast_problem.b+elast_problem.b_frac);
R=elast_problem.R;

f=elast_problem.B_K_plus*b-elast_problem.c;
e=R'*b;

lambda_ImGt=elast_problem.G'*((elast_problem.GGt)\(e));

c_ker=-lambda_ImGt;
c_ker(elast_problem.idx_no_bounds)=0;

%%
update_temp_struct.lambda_ImGt=lambda_ImGt;

update_G=elast_problem.update_G;

rel=SMALSE_params.rel;
rho0=SMALSE_params.rho0;
betarho=SMALSE_params.betarho;
Gama=SMALSE_params.Gama;
M_start=SMALSE_params.M_start;
tol_to_update=SMALSE_params.tol_to_update;
maxiter_cg=SMALSE_params.maxiter_cg;
type=SMALSE_params.type;
print=SMALSE_params.print;
update_G_empty=[];
if isempty(elast_problem.lambda_ker)
    [lambda_ker,mi,update_temp_struct,~,ncg] = smalse.SMALSE_update_new...
        (elast_problem.F,f-elast_problem.F*lambda_ImGt,elast_problem.G,elast_problem.Q,c_ker,update_G_empty,[],[],update_temp_struct,elast_problem.idx_no_bounds,...
        rel,tol_to_update,rho0,betarho,Gama,M_start,maxiter_cg,type,print,elast_problem.norms_mat);
else
    [lambda_ker,mi,update_temp_struct,~,ncg] = smalse.SMALSE_update_new...
        (elast_problem.F,f-elast_problem.F*lambda_ImGt,elast_problem.G,elast_problem.Q,c_ker,update_G_empty,elast_problem.lambda_ker,elast_problem.mi,update_temp_struct,elast_problem.idx_no_bounds,...
        rel,tol_to_update,rho0,betarho,Gama,M_start,maxiter_cg,type,print,elast_problem.norms_mat);
end
elast_problem.lambda_ker=lambda_ker;
elast_problem.mi=mi;
[x_elast] =update_G(lambda_ker,b,update_temp_struct);

[D] = a_ela.construct_apertures(elast_problem.B_i*(x_elast+elast_problem.NODES),elast_problem.fracture_matrice);
elast_problem.x_elast=x_elast;
end

