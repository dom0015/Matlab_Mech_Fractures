%% load pre-constructed problem
Input_Data=load('data_two_fractures.mat');

%% linear elasticity with TFETI (B_i) and contact in fractures (B_e)
ELEMS_local=Input_Data.ELEMS_local;
NODES_local=Input_Data.NODES_local;

K=Input_Data.K;
K_plus=Input_Data.K_plus;
R=Input_Data.R;

k=Input_Data.b;

B_e=Input_Data.B_e;
B_i=Input_Data.B_i;

c_e=Input_Data.c_e;
c_i=Input_Data.c_i;

%% dualization
c=[c_e;c_i];
B_ei=[B_e;-B_i];

F=B_ei*K_plus*B_ei';
G=R'*B_ei';
f=B_ei*K_plus*k-c;
e=R'*k;

%% homogenization
lambda_ImGt=G'*((G*G')\e);

%% orthonormalize eq. constraint
L_GG=chol(G*G')';
%L_GG=speye(size(G*G'));
G_fnc=@(x)L_GG\(G*x);
Gt_fnc=@(x)G'*(L_GG'\x);

%% final system min(0.5xAx-xb) over Gx=0 and x>=lb
A=F;
b=f-F*lambda_ImGt;

scale=eigs(A,1);
A=A/scale;
b=b/scale;

%% Project onto natural coarse space
P_fnc=@(x)x-Gt_fnc(G_fnc(x));
P_fnc=@(x)x;
A_fnc=@(x)P_fnc(A*P_fnc(x));
b=P_fnc(b);

lb=-lambda_ImGt;
lb(1:length(c_e))=-inf;


%% solver

mprgpctx.lb = lb;
mprgpctx.abarmult = 2;
mprgpctx.propConst = 1;
mprgpctx.settol = 10*eps;
mprgpctx.infeastol = 0;
type='M';
tic;
[x,flg,k,mu,all_errors,all_iters] = smalxe_fnc(A_fnc,b,G_fnc,Gt_fnc,zeros(size(G,1),1),zeros(size(A,1),1),zeros(size(G,1),1),1,1.1,1,1.5,1e-1,1e-6,1000,@mprgp_fnc,type,10000,mprgpctx);
toc
plot_convergence(all_iters,all_errors)
sum(all_iters)

tic;
[x,flg,k,mu,all_errors,all_iters] = smalxe_fnc(A_fnc,b,G_fnc,Gt_fnc,zeros(size(G,1),1),x,mu,1,1.1,1,1.5,1e-1,1e-12,1000,@mprgp_fnc,type,10000,mprgpctx);
toc
plot_convergence(all_iters,all_errors)
sum(all_iters)
