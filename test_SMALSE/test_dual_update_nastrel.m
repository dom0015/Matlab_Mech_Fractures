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
relgeom=1e-8;
relupd=1e-12;
rho0=2;
betarho=1.5;
Gama = 1;
M_start=0.1;
tol_to_update=1e-8;
maxiter_cg = 1000;
type='m';
prnt=true;

tic
[lambda_ker,update_temp_struct] = SMALSE_update(F,d-F*lambda_ImGt,G,c_ker,[],update_G,update_temp_struct,idx_no_bounds,rel,relgeom,relupd,tol_to_update,rho0,betarho,Gama,M_start,maxiter_cg,type,prnt);
toc

[F,ff,G,c_ker,update_temp_struct,geom_res,x_elast] =update_G(lambda_ker,update_temp_struct);

plot_func(x_elast,fracture_matrice);

% rel=1.0e-13;
% relgeom=1e-9;
% 
% tic
% [lambda_ker,update_temp_struct] = SMALSE_update(F,ff,G,c_ker,lambda_ker,update_G,update_temp_struct,idx_no_bounds,rel,relgeom,relupd,tol_to_update,rho0,betarho,Gama,M_start,maxiter_cg,type,prnt);
% toc
