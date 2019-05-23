
mat_omega = 1e-15; % material - matrice
mat_fracture = 1e-9;
d = exp([-6 -8.5]);
alfa_intersections = 1e-14;
p_Dirichlet=1e6; % pressure on Dir. b. c.

[u0,freeNode,B,b,elem,node,A,fractures_positions,lengths] = wrap_cg(mat_omega,mat_fracture,d,alfa_intersections,p_Dirichlet);

Bf=B(freeNode,freeNode);
bf=b(freeNode);

uf_=Bf\bf;

L=chol(Bf,'lower');
%L=ichol(Bf);
[uf,flag,relres,iter,resvec]=pcg(Bf,bf,1e-16,length(bf),L,L');
%[uf,flag,relres,iter,resvec]=gmres(Bf,bf,[],1e-9,300);

figure

plot(resvec)
set(gca,'YScale','log')
norm(uf_-uf)/norm(uf_)

u0(freeNode)=uf_;
wrap_cg_figure(A,node,elem,u0,fractures_positions,p_Dirichlet,lengths);