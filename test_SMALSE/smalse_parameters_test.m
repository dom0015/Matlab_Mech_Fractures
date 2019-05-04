% d_eps=1e-8;
% d_min=0;
% max_it=1000;
% rho_step=1.2;
% rho_step_limits=10;
% stagnate_count=20;
% alpha=0.75;
% 
% [u,i,d_vec,un_vec,geom_vec,rho] = penalta_FETI_noswitch(A,B_e,b,c_e,B_iupdate,d_min,max_it,d_eps,rho_step,rho_step_limits,stagnate_count,alpha);
% 
% B_i=B_iupdate(u);


B=[B_e;-B_i];
c=[c_e;c_i];
F=B*A_plus*B';
d=B*A_plus*b;
G=R'*B';
e=R'*b;

lambda_ImGt=G'*((G*G')\(e));
idx_no_bounds=1:length(c_e);
idx_bounds=length(c_e)+1:length(c);
c_ker=-lambda_ImGt;
c_ker(idx_no_bounds)=0;


update_temp_struct.lambda_ImGt=lambda_ImGt;
update_temp_struct.B_i=B_i;
update_temp_struct.Be_Ap_Be=B_e*A_plus*B_e';
update_temp_struct.Be_Ap=B_e*A_plus;
update_temp_struct.Ap_b=A_plus*b;



update_G=@(x,data)Update_all_dual(x,idx_no_bounds,idx_bounds,c,B_e,A_plus,b,R,B_iupdate,data);





rel=1.0e-12;
rho0=10;
betarho=2;
Gama = 1;
M_start=0.1;
tol_to_update=1e-7;
maxiter_cg = 10000;
type='M';


tic
[lambda_ker,update_temp_struct] = SMALSE_update(F,d-F*lambda_ImGt,G,c_ker,update_G,update_temp_struct,idx_no_bounds,rel,tol_to_update,rho0,betarho,Gama,M_start,maxiter_cg,type,true);
toc

[~,~,~,~,~,~,x_elast] =update_G(lambda_ker,update_temp_struct);


plot_func(x_elast,fracture_matrice);

