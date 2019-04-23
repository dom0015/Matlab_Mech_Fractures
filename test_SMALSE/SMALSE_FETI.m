%% PARAMETERS -------------------------------------------------------------
Nxy=101;
L1=1; L2=1;
sumbdomains_FETI=ceil(Nxy/15)^2;
frac_start_end={[0.8 0.1], [0.8 0.9]
     [0.1 0.2], [0.9 0.2]
     [0.1 0.1], [0.9 0.9]
     [0.1 0.5], [0.9 0.5]};


% frac_start_end={[0.5 0.1], [0.5 0.9]
%      [0.1 0.4], [0.9 0.4]
%      [0.1 0.8], [0.9 0.8]
%      [0.1 0.2], [0.9 0.2]
%      [0.1 0.1], [0.9 0.9]};
FETI_test_assembly



%A,b=cat(1,b{:});
% A_plus = blkdiag(A_pinv{:});
% R = blkdiag(A_null{:});
% B_e=[B_dirichlet;B_FETI];
% B_i = contact_inequalities(0*b_full,POINTS,fracture_matrice,node_map_on,size(A,1));
% c_e=[B_dirichlet_rhs;zeros(size(B_FETI,1),1)];
% c_i=zeros(size(B_i,1),1);
% 
% B_iupdate =@(x) contact_inequalities(x,POINTS,fracture_matrice,node_map_on,size(A,1));

%B_e= my_orth(B_e,1e-16);
B=[B_e;-B_i];
c=[c_e;c_i];
F=B*A_plus*B';
d=B*A_plus*b;
G=R'*B';
e=R'*b;

lambda_0=[x_full(size(A,1)+1:end);c_i];
lambda_0=lambda_0*0;

lambda_ImGt=G'*((G*G')\(e));
idx_no_bounds=1:length(c_e);
idx_bounds=length(c_e)+1:length(c);
c_ker=-lambda_ImGt;
c_ker(idx_no_bounds)=0;


%%
rel=1.0e-8;
rho0=1.2;
betarho=1.01;
Gama = 1;
maxiter_cg = 10000;
%%

[lambda_ker] = SMALSE(F,d-F*lambda_ImGt,G,c_ker,idx_no_bounds,rel,rho0,betarho,Gama,maxiter_cg);
lambda=lambda_ker+lambda_ImGt;
x_elast=A_plus*(b-B'*lambda);

%alpha=(B_e*R)\(B_e*x_elast);

idx_active_equality=false(length(c),1);
idx_active_equality(idx_no_bounds)=true;
idx_active_equality(idx_bounds)=lambda(idx_bounds)>0;
B_pruh=B(idx_active_equality,:);
c_pruh=c(idx_active_equality);
alpha=pinv(full(R'*B_pruh'*B_pruh*R))*R'*B_pruh'*(c_pruh-B_pruh*x_elast);
x_elast=x_elast+R*alpha;








%%%%
% for i=1:1
% B_i=B_iupdate(x_elast);
% B=[B_e;-B_i];
% c=[c_e;c_i];
% F=B*A_plus*B';
% d=B*A_plus*b;
% G=R'*B';
% e=R'*b;
% 
% lambda_0=[x_full(size(A,1)+1:end);c_i];
% lambda_0=lambda_0*0;
% 
% lambda_ImGt=G'*((G*G')\(e));
% idx_no_bounds=1:length(c_e);
% idx_bounds=length(c_e)+1:length(c);
% c_ker=-lambda_ImGt;
% c_ker(idx_no_bounds)=0;
% 
% 
% %%
% rel=1.0e-8;
% rho0=1;
% betarho=1.01;
% Gama = 1;
% maxiter_cg = 10000;
% %%
% 
% [lambda_ker] = SMALSE(F,d-F*lambda_ImGt,G,c_ker,idx_no_bounds,rel,rho0,betarho,Gama,maxiter_cg);
% lambda=lambda_ker+lambda_ImGt;
% x_elast=A_plus*(b-B'*lambda);
% 
% %alpha=(B_e*R)\(B_e*x_elast);
% 
% idx_active_equality=false(length(c),1);
% idx_active_equality(idx_no_bounds)=true;
% idx_active_equality(idx_bounds)=lambda(idx_bounds)>0;
% B_pruh=B(idx_active_equality,:);
% c_pruh=c(idx_active_equality);
% alpha=(R'*(B_pruh'*B_pruh)*R)\(R'*B_pruh'*(c_pruh-B_pruh*x_elast));
% 
% [alpha_,flag,relres,iter,resvec]=gmres((R'*(B_pruh'*B_pruh)*R),(R'*B_pruh'*(c_pruh-B_pruh*x_elast)),[],1e-12,100);
% x_elast=x_elast+R*alpha_;
% end


plot_func(x_elast*1000,fracture_matrice);

% B_i_onidx=lambda(idx_bounds)>0;
% B_i_on=B_i(B_i_onidx,:);
% B_lambda=[B_e;B_i_on];
% M=[A B_lambda';B_lambda sparse(size(B_lambda,1),size(B_lambda,1))];
% f=[b;c_e;c_i(B_i_onidx)];
% x_lambda=M\f;
% 
% 
% plot_func(x_lambda,fracture_matrice);



