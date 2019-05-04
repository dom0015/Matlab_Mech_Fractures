load('test_stability.mat')
idx_no_bounds=1:length(c_e);
idx_bounds=length(c_e)+1:length(c_i);

%% parametry SMALSE--------------------------------------------------------
rel=1.0e-14;
rho0=1;
betarho=2;
Gama = 1;
maxiter_cg = 10000;
M_start=0.1;
type='M';

%perturbace
perturbace=1e-12; % velikost perturbace posunuti -> vede ke drobne zmene normalv trhlinach ~ B_i

%% puvodni B_i ->B_i je slozene ze dvou stejnych matic za sebou ~ normaly z leva a z prava jsou stejne
B_i=B_iupdate(b*0); % normály oproti nulovemu posunu

%sestaveni matic pro dualni formulaci
B=[B_e;-B_i];
c=[c_e;c_i];
F=B*A_plus*B';
G=R'*B';
d=B*A_plus*b;
e=R'*b;
%partikularni reseni lambda_img
lambda_ImGt=G'*((G*G')\(e));
c_ker=-lambda_ImGt;
c_ker(idx_no_bounds)=0;

%SMALSE reseni
[lambda_ker] = SMALSE(F,(d-F*lambda_ImGt),G,c_ker,idx_no_bounds,rel,rho0,betarho,Gama,M_start,maxiter_cg,type);

%rekonstrukce posunuti
lambda=lambda_ker+lambda_ImGt;
x_elast=A_plus*(b-B'*lambda);

idx_active_equality=false(length(c),1);
idx_active_equality(idx_no_bounds)=true;
idx_active_equality(idx_bounds)=lambda(idx_bounds)>0;
B_pruh=B(idx_active_equality,:);
c_pruh=c(idx_active_equality);
B_R=B_pruh*R;
alpha=(B_R'*B_R)\(B_R'*(c_pruh-B_pruh*x_elast)); %sestaveni alpha jako indexy doplnku x v kernelu A

x_elast=x_elast+R*alpha;

%vykresleni
plot_func(x_elast,fracture_matrice);

%% perturbovane B_i drobne rozdily mezi normálami zleva a zprava ----------

x_tmp=0*b+(rand(length(A),1)-0.5)*2*perturbace; % perturbace ~ miniaturni posuny bodu
B_i_perturb=B_iupdate(x_tmp); % normaly oproti miniaturnim posunum

%sestaveni matic pro dualni formulaci
B=[B_e;-B_i_perturb];
c=[c_e;c_i];
F=B*A_plus*B';
G=R'*B';
d=B*A_plus*b;
e=R'*b;
%partikularni reseni lambda_img
lambda_ImGt=G'*((G*G')\(e));
c_ker=-lambda_ImGt;

%SMALSE reseni
[lambda_ker] = SMALSE(F,(d-F*lambda_ImGt),G,c_ker,idx_no_bounds,rel,rho0,betarho,Gama,M_start,maxiter_cg,type);

%rekonstrukce posunuti
lambda=lambda_ker+lambda_ImGt;
x_elast=A_plus*(b-B'*lambda);

idx_active_equality=false(length(c),1);
idx_active_equality(idx_no_bounds)=true;
idx_active_equality(idx_bounds)=lambda(idx_bounds)>0;
B_pruh=B(idx_active_equality,:);
c_pruh=c(idx_active_equality);
B_R=B_pruh*R;
alpha=(B_R'*B_R)\(B_R'*(c_pruh-B_pruh*x_elast)); %sestaveni alpha jako indexy doplnku x v kernelu A

x_elast=x_elast+R*alpha;

%vykresleni
plot_func(x_elast,fracture_matrice);
