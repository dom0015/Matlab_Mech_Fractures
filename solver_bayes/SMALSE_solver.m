function [D,problem_setting] = SMALSE_solver(problem_setting,SMALSE_params)
%SMALSE_SOLVER Summary of this function goes here
%   Detailed explanation goes here

A_plus=problem_setting.A_plus*problem_setting.mat_scale;
B_e=problem_setting.B_e;
b=(problem_setting.b+problem_setting.b_frac)/problem_setting.mat_scale;
R=problem_setting.R;
c_e=problem_setting.c_e;
c_i=problem_setting.c_i;
B_iupdate=problem_setting.B_iupdate;
fracture_matrice=problem_setting.fracture_matrice;
plot_func2=problem_setting.plot_func2;
if isempty(problem_setting.B_i)
    B_i=B_iupdate(0*b);
else
    B_i=problem_setting.B_i;
end
B=[B_e;-B_i];
c=[c_e;c_i*ones(size(B_i,1),1)];

F=B*A_plus*B';
G=R'*B';
d=B*A_plus*b;

e=R'*b;

lambda_ImGt=G'*((G*G')\(e));
idx_no_bounds=1:length(c_e);
idx_bounds=length(c_e)+1:length(c);
c_ker=-lambda_ImGt;
c_ker(idx_no_bounds)=0;

%%
update_temp_struct.lambda_ImGt=lambda_ImGt;
update_temp_struct.B_i=B_i;
update_temp_struct.Be_Ap_Be=B_e*A_plus*B_e';
update_temp_struct.Be_Ap=B_e*A_plus;
update_temp_struct.Ap_b=A_plus*b;

update_G=@(x,data)Update_all_dual(x,idx_no_bounds,idx_bounds,c,B_e,A_plus,b,R,B_iupdate,data);

rel=SMALSE_params.rel;
rho0=SMALSE_params.rho0;
betarho=SMALSE_params.betarho;
Gama=SMALSE_params.Gama;
M_start=SMALSE_params.M_start;
tol_to_update=SMALSE_params.tol_to_update;
maxiter_cg=SMALSE_params.maxiter_cg;
type=SMALSE_params.type;
print=SMALSE_params.print;
if isempty(problem_setting.lambda_ker)
[lambda_ker,mi,update_temp_struct] = SMALSE_update...
    (F,d-F*lambda_ImGt,G,c_ker,update_G,[],[],update_temp_struct,idx_no_bounds,...
    rel,tol_to_update,rho0,betarho,Gama,M_start,maxiter_cg,type,print);
else
    [lambda_ker,mi,update_temp_struct] = SMALSE_update...
    (F,d-F*lambda_ImGt,G,c_ker,update_G,problem_setting.lambda_ker,problem_setting.mi,update_temp_struct,idx_no_bounds,...
    rel,tol_to_update,rho0,betarho,Gama,M_start,maxiter_cg,type,print);
end
problem_setting.B_i=update_temp_struct.B_i;
problem_setting.lambda_ker=lambda_ker;
problem_setting.mi=mi;
[~,~,~,~,~,~,x_elast] =update_G(lambda_ker,update_temp_struct);

if print
plot_func2(x_elast,fracture_matrice);
end
[D] = construct_apertures(update_temp_struct.B_i*x_elast,fracture_matrice);
problem_setting.x_elast=x_elast;
end

