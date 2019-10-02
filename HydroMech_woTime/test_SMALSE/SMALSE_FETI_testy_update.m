%% PARAMETERS -------------------------------------------------------------
Nxy=101;
L1=1; L2=1;
sumbdomains_FETI=ceil(Nxy/15)^2;


% mat_const=1;
% frac_start_end={[0.1 0.5], [0.9 0.5]};
% frac_press={@(x)(-0.5+sin(5*pi*x)),@(x)0*x};



% mat_const=10000000;
% frac_start_end={[0.1 0.4], [0.9 0.4]
%      [0.3 0.5], [0.6 0.5]
%      [0.1 0.6], [0.9 0.6]};
% frac_press={@(x)mat_const/10+0*x,@(x)mat_const/10+0*x
%     @(x)2*mat_const/10+0*x,@(x)2*mat_const/10+0*x
%     @(x)mat_const/10+0*x,@(x)mat_const/10+0*x};
% FETI_test_assembly

% mat_const=10;
% frac_start_end={[0.1 0.4], [0.9 0.4]
%      [0.3 0.5], [0.6 0.5]
%      [0.1 0.6], [0.9 0.6]};
% frac_press={@(x)1+0*x,@(x)1+0*x
%     @(x)2+0*x,@(x)2+0*x
%     @(x)1+0*x,@(x)1+0*x};
%FETI_test_assembly


mat_const=1;
frac_start_end={[0.1 0.4], [0.9 0.4]
     [0.3 0.5], [0.6 0.5]
     [0.1 0.6], [0.9 0.6]
     [0.5 0.1], [0.5 0.9]};
frac_press={@(x)mat_const/5+0*x,@(x)mat_const/5+0*x
    @(x)2*mat_const/5+0*x,@(x)2*mat_const/5+0*x
    @(x)mat_const/5+0*x,@(x)mat_const/5+0*x
    @(x)mat_const/5+0*x,@(x)mat_const/5+0*x};

FETI_test_assembly

global lambda_ImGt B_i
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


%%
rel=1.0e-12;
rho0=10;
betarho=1.5;
Gama = 1;
maxiter_cg = 10000;
%%

update_G=@(x)Update_B_dual(x,idx_no_bounds,idx_bounds,c,B_e,A_plus,b,R,B_iupdate);

[lambda_ker] = SMALSE_dual_update(F,d-F*lambda_ImGt,G,c_ker,idx_no_bounds,rel,rho0,betarho,Gama,maxiter_cg,update_G);

lambda=lambda_ker+lambda_ImGt;
B=[B_e;-B_i];
x_elast=A_plus*(b-B'*lambda);

% sestaveni alpha jako indexy doplnku x v kernelu A
idx_active_equality=false(length(c),1);
idx_active_equality(idx_no_bounds)=true;
idx_active_equality(idx_bounds)=lambda(idx_bounds)>0;
B_pruh=B(idx_active_equality,:);
c_pruh=c(idx_active_equality);
B_R=B_pruh*R;
alpha=(B_R'*B_R)\(B_R'*(c_pruh-B_pruh*x_elast));

x_elast=x_elast+R*alpha;


plot_func(x_elast,fracture_matrice);