function [D,elast_problem,x_elast,ncg] = SMALSE_solver(elast_problem,SMALSE_params)
%SMALSE_SOLVER Summary of this function goes here
%   Detailed explanation goes here
K_plus=elast_problem.A_plus;%*problem_setting.mat_scale;
B_e=elast_problem.B_e;
b=(elast_problem.b+elast_problem.b_frac);%/problem_setting.mat_scale;
R=elast_problem.R;
c_e=elast_problem.c_e;
B_iupdate=elast_problem.B_iupdate;
fracture_matrice=elast_problem.fracture_matrice;
plot_func2=elast_problem.plot_func2;
if isempty(elast_problem.B_i)
    B_i=B_iupdate(0*b);
else
    B_i=elast_problem.B_i;
end
NODES = cell2mat(elast_problem.sub_nodes)';
NODES = NODES(:);
c_i=B_i*NODES;
if ~isempty(elast_problem.c_i)
    c_i=c_i-elast_problem.c_i;
end
c=[c_e;c_i];
B=[B_e;-B_i];

F=B*K_plus*B';
G=R'*B';
f=B*K_plus*b-c;
e=R'*b;

lambda_ImGt=G'*((G*G')\(e));

idx_no_bounds=1:length(c_e);
idx_bounds=length(c_e)+1:length(c);
c_ker=-lambda_ImGt;
c_ker(idx_no_bounds)=0;

%%
update_temp_struct.lambda_ImGt=lambda_ImGt;
update_temp_struct.B_i=B_i;
update_temp_struct.Be_Ap_Be=B_e*K_plus*B_e';
update_temp_struct.Be_Ap=B_e*K_plus;
update_temp_struct.Ap_b=K_plus*b;

update_G=@(x,data)smalse.Update_all_dual(x,idx_no_bounds,idx_bounds,c,B_e,K_plus,b,R,B_iupdate,data);

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
    [lambda_ker,mi,update_temp_struct,~,ncg] = smalse.SMALSE_update...
        (F,f-F*lambda_ImGt,G,c_ker,update_G_empty,[],[],update_temp_struct,idx_no_bounds,...
        rel,tol_to_update,rho0,betarho,Gama,M_start,maxiter_cg,type,print);
else
    [lambda_ker,mi,update_temp_struct,~,ncg] = smalse.SMALSE_update...
        (F,f-F*lambda_ImGt,G,c_ker,update_G_empty,elast_problem.lambda_ker,elast_problem.mi,update_temp_struct,idx_no_bounds,...
        rel,tol_to_update,rho0,betarho,Gama,M_start,maxiter_cg,type,print);
end
elast_problem.B_i=update_temp_struct.B_i;
elast_problem.lambda_ker=lambda_ker;
elast_problem.mi=mi;
[~,~,~,~,~,~,x_elast] =update_G(lambda_ker,update_temp_struct);

[D] = a_ela.construct_apertures(update_temp_struct.B_i*(x_elast+NODES),fracture_matrice);
elast_problem.x_elast=x_elast;
if print
    plot_func2(x_elast,fracture_matrice);
    figure(103); d_viz = cell2mat(D); d_viz(1)=d_viz(1)*2; d_viz(end)=d_viz(end)*2; plot(d_viz); title('All apertures from SMALSE'); hold on
end

end

