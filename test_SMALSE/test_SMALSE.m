% SMALSE David
% MPRGP

%% inputs (from elesticity with fracture)
n_node=101; % 5*k+1;
[F,b0,G,geom] = SMALSE_precomp(n_node);

[m,n]=size(G);
F=[F zeros(n,m); zeros(m,n+m)];
b0=[b0; zeros(m,1)];
c=0*b0;
idx_no_bounds=1:n;
G=[G eye(m)];

update_G=@(x)[geom.recalculate_B(x(idx_no_bounds)) eye(m)];
%%

rel=1.0e-10;
rho0=1;
betarho=1.1;
Gama = 1;
maxiter_cg = 100000;


[u] = SMALSE(F,b0,G,c,idx_no_bounds,rel,rho0,betarho,Gama,maxiter_cg,update_G);
G_new=geom.recalculate_B(u(1:n));
G=[G_new eye(m)];
figure
semilogy(abs((u(n+1:end)+G_new*u(1:n))))


geom.plot(u(1:n));
