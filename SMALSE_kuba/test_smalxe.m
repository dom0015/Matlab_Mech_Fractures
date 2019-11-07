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

scale=eigs(A,1);
A=A/scale;
b=b/scale;

B=G;
lb=-lambda_ImGt;
lb(1:length(c_e))=-inf;


%% solver

mprgpctx.lb = lb;
mprgpctx.abarmult = 2;
mprgpctx.propConst = 1;
mprgpctx.settol = 10*eps;
mprgpctx.infeastol = 0;
tic;
[x,flg,k] = smalxe(A,b,B,zeros(size(B,1),1),zeros(size(A,1),1),10,2,10,2,1e-1,1e-10,100,@mprgp,1000,mprgpctx);
toc

