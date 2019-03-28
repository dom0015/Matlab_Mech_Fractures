% SMALSE David
% MPRGP

%% inputs (from elesticity with fracture)
n_node=101; % 5*k+1;
[F,b0,G,geom] = SMALSE_precomp(n_node);

u_p = penalta(F,b0,geom.recalculate_B,0,100,1e-7,2);
G=geom.recalculate_B(u_p);
slack_=-G*u_p;
slack_(slack_<0)=0;

geom.plot(u_p);

[m,n]=size(G);
F=[F zeros(n,m); zeros(m,n+m)];
b0=[b0; zeros(m,1)];
c=0*b0;
idx_no_bounds=1:n;
G=[G eye(m)];
Q=G'*G;
%%

rel=1.0e-12;  
rho0=1; 
betarho=1.1;
Gama = 1;
maxiter_cg = 100000;

[u] = SMALSE(F,b0,G,Q,c,idx_no_bounds,rel,rho0,betarho,Gama,maxiter_cg);

geom.plot(u(1:n));
