% SMALSE David
% MPRGP

%% inputs (from elesticity with fracture)
n_node=41; % 5*k+1;
[A,b0,B,geom] = SMALSE_precomp(n_node);
u_p = penalta(A,b0,geom.recalculate_B,0,100,1e-7,2);
geom.plot(u_p);
B=geom.recalculate_B(u_p);


[m,n]=size(B);
clear opts
%opts.type = 'nofill';
opts.michol = 'on';
 %L=ichol(A+1000*((B')*B),opts);
L=chol(A,'lower');
Ct_inv=@(x)L\x;
C_inv=@(x)(L')\x;

A_mat=[speye(n) zeros(n,m); zeros(m,n+m)];
BtB=B'*B;

b0=[Ct_inv(b0); zeros(m,1)];
c=0*b0;
idx_no_bounds=1:n;

F=@(x)A_mat*x;
%F=@(x)[Ct_inv(A*C_inv(x(1:n)));zeros(m,1)];
update_GQ=@(x)update_precond_G(x,geom.recalculate_B,C_inv,Ct_inv,n);

%%

rel=1.0e-10;
rho0=200;
betarho=1.8;
Gama = 1;
maxiter_cg = 100000;


%[u] = SMALSE(F,b0,G,c,idx_no_bounds,rel,rho0,betarho,Gama,maxiter_cg,update_G);

[u] = SMALSE_matfun(F,n,m,b0,c,idx_no_bounds,rel,rho0,betarho,Gama,maxiter_cg,update_GQ);
u(1:n)=C_inv(u(1:n));
G_new=geom.recalculate_B(u(1:n));
G=[G_new eye(m)];
figure
semilogy(abs((u(n+1:end)+G_new*u(1:n))))


geom.plot(u(1:n));
