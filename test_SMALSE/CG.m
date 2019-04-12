function [ u, res ] = CG( A,P,b,u,epsilon,maxit )
% Metoda sdruzenych gradientu (Conjugate gradient method),
% pouziva pouze dva skalarni souciny.
if nargin<2
    b=ones(length(A),1);
end
if nargin<3
    u=zeros(length(b),1);
end
if nargin<4
    epsilon=1e-6;
end
if nargin<5
    maxit=100;
end

r=b-A*P(u);
d=r;
rho0=r'*r;
rho_=rho0;
res=[];
for i=1:maxit
    w=A*P(d);
    alpha=rho_/(w'*d);
    u=u+alpha*d;
    r=r-alpha*w;
    rho=r'*r;
    res=[res rho/rho0];
    if rho<=epsilon^2*rho0
        break
    end
    beta=rho/rho_;
    d=r+beta*d;
    rho_=rho;
end

end

