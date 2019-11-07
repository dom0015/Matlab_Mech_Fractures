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


%% final system min(0.5xAx-xb) over Gx=0 and x>=lb
A=F;
b=f-F*lambda_ImGt;
B=G;
lb=-lambda_ImGt;
lb(1:length(c_e))=-inf;

%% solver
rel=1e-10;
rho0=1;
betarho=2;
Gama=2;
M0=1;
maxiter_cg=10000;
type='m';
printres=1;
tic;
[u,mi,nout,ncg] = old_smalse(A,b,B,lb,[],[],rel,rho0,betarho,Gama,M0,maxiter_cg,type,printres);
toc