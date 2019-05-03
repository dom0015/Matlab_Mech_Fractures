
Lt_inv=inv(chol(B_e*B_e')');
B_e_orth=B_e;
c_e_orth=c_e;
B=[B_e_orth;-B_i];
c=[c_e_orth;c_i];
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
rel=1.0e-14;
rho0=1;
betarho=2;
Gama = 1;
maxiter_cg = 10000;
M_start=0.1;
type='M';

[lambda_ker] = SMALSE_nastrel(F,(d-F*lambda_ImGt),G,c_ker,[],idx_no_bounds,rel,rho0,betarho,Gama,M_start,maxiter_cg,type);

%% perturbovane

B_e_orth=B_e;
c_e_orth=c_e;

x_tmp=zeros(length(A),1)+(rand(length(A),1)-0.5)*1e-10;
B_i_perturb=B_iupdate(x_tmp);

B=[B_e_orth;-B_i_perturb];
Lt_inv=inv(chol(B*B'+speye(size(B,1))*1e-10)');
B_orth=Lt_inv*B;


c=[c_e;c_i];
c=Lt_inv*c;
F=B_orth*A_plus*B_orth';
G=R'*B_orth';
d=B_orth*A_plus*b;

e=R'*b;

lambda_ImGt=G'*((G*G')\(e));
idx_no_bounds=1:length(c_e);
idx_bounds=length(c_e)+1:length(c);
c_ker=-lambda_ImGt;
c_ker(idx_no_bounds)=0;



%%
rel=1.0e-14;
rho0=1;
betarho=2;
Gama = 1;
maxiter_cg = 10000;
M_start=0.1;
type='M';

[lambda_ker] = SMALSE_nastrel(F,(d-F*lambda_ImGt),G,c_ker,[],idx_no_bounds,rel,rho0,betarho,Gama,M_start,maxiter_cg,type);
