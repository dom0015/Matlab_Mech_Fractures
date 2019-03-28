% SMALSE David
% MPRGP

%% inputs (from elesticity with fracture)
n_node=41; % 5*k+1;
[A,b0,B,geom] = SMALSE_precomp(n_node);

u_p = penalta(A,b0,geom.recalculate_B,0,100,1e-7,2);
B=geom.recalculate_B(u_p);

F=B*(A\(B'));
u0=

%R=null(A);


rel=1.0e-12;
rho0=1;
betarho=1.1;
Gama = 1;
maxiter_cg = 100000;





[u] = SMALSE(F,b0,G,c,idx_no_bounds,rel,rho0,betarho,Gama,maxiter_cg);
G_new=geom.recalculate_B(u(1:n));
G=[G_new eye(m)];
figure
semilogy(abs((u(n+1:end)+G_new*u(1:n))))


geom.plot(u(1:n));
