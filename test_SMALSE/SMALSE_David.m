% SMALSE David
% MPRGP

%% inputs (from elesticity with fracture)
n_node=81; % 5*k+1;
[F,b0,G,geom] = SMALSE_precomp(n_node);

u = penalta(F,b0,geom.recalculate_B,0,100,1e-4,2);
G=geom.recalculate_B(u);

[m,n]=size(G);
F=[F zeros(n,m); zeros(m,n+m)];
b0=[b0; zeros(m,1)];
c=0*b0;
niq=1:n;
G=[G eye(m)];
Q=G'*G;
%%

lFl=normest(F);


rel=1.0e-4;  rho0=10; betarho=2;
epsr = rel*norm(b0); Gama = 1;
ncg = 0; ne = 0; np = 0; nout = 0; ncg_max = 10000;
VLag=[]; VLagIn=[]; Vnarus=[]; Vgp=[]; VM=[]; VCrit=[]; Increment=0;

c(niq) = -Inf*ones(length(niq),1);
%u = max(B*f/2,c);
u=max(zeros(size(F,1),1),c);
J = (u > c);
rho = rho0* lFl;
Gu=G*u;
mi=zeros(size(G,1),1);  mi=mi+rho*Gu;
b=b0-G'*mi;   b_old=b;
% A=P*F*P+rho*Q;
A=F+rho*Q;
lAl=normest(A);
alfa=1/lAl;         % alpha_bar...can be computed more precisously  by Reigh.quotient
 M=10*lAl; M_old = M;
g = A*u - b;

Lag=0.5*u'*A*u-b'*u;
gp=ones(length(u),1);
ny0=norm(b0);

while (1)
    %geom.plot(u(1:n));
    if (norm(gp)<epsr && norm(Gu)<epsr) break; end;
    
    %     lAl=norm(full(A));
    alfafix=1/lAl;
    
    % Gradient splitting of g=-r, gf=gradient free(fi), gr=gradient free
    % reduced (fired), gc=gradient chopped or cut (beta), gp=gf+gc=projected grad.
    % previously used gr = min(lAl*J.*(u-c), gf);
    g_old=g;
    g = A*u - b; gf = J.*g; gc = min((~J).*g,0); gr = min(J.*(u-c)/alfafix, gf); gp = gf + gc;
    p=gf;      Gu_old=Gu;       Gu=G*u;
    
    while (1)
        crit=min(M*norm(Gu),ny0);
        %fprintf('%f < 1 nebo %f<1 a %f < 1\n',norm(gp)/ crit,norm(gp)/epsr, norm(Gu)/epsr)
        if norm(gp)<crit || (norm(gp)<epsr && norm(Gu)<epsr) break; end;
        VCrit=[VCrit crit];
        
        if gc'*gc <= Gama^2*gr'*gf
            % Proportional iteration. Trial conjugate gradient step.
            Ap=A*p; rtp=g'*p; pAp=p'*Ap; acg=rtp/pAp;  yy=u-acg*p; ncg=ncg+1;
            if all(yy>=c)
                %Conjugate gradient step
                u=yy; g=g-acg*Ap;
                gf = J.*g; gc = min((~J).*g,0); gr = min(lAl*J.*(u-c), gf); gp = gf + gc;
                beta=gf'*Ap/pAp;
                p=gf-beta*p;
            else
                %partial conjugate gradient step
                if max((J.*p)./(J.*(u-c)+(~J)))~=0
                    a=1/max((J.*p)./(J.*(u-c)+(~J)));
                else
                    a=0;
                end
                u=max(u-a*p,c); J=(u>c); g=g-a*Ap;
                gf = J.*g; gr = min(lAl*J.*(u-c), gf);
                %expansion step
                %                 alfa = (gr'*g)/(gr'*A*gr);
                alfa = 1/lAl;
                u=max(u-alfa*gr,c); J=(u>c); g=A*u-b;
                gf = J.*g; gc = min((~J).*g,0); gr = min(lAl*J.*(u-c), gf); gp = gf + gc;
                p=gf; ne=ne+1;
            end
        else
            %Proportioning step
            Ap=A*gc; acg=(gc'*g)/(gc'*Ap);
            u=u-acg*gc; J=(u>c); g=g-acg*Ap;
            gf = J.*g; gc = min((~J).*g,0); gr = min(lAl*J.*(u-c), gf); gp = gf + gc;
            p=gf; np=np+1;
        end
        
        Gu=G*u;  VLagIn=[VLagIn 0.5*u'*(g-b)]; Vnarus=[Vnarus norm(Gu)]; Vgp=[Vgp norm(gp)];
        if (ncg>=ncg_max) break; end;
    end
    
    disp(sprintf('number of CG iterations = %d', ncg));
    mi_old=mi;        mi=mi+rho*Gu;
    VLag=[VLag Lag];  VM=[VM M];
    Lag_old=Lag;      Lag=0.5*u'*A*u-b'*u;
    
    Increment=0.5*rho*crit^2;
    
    if Lag - Lag_old < Increment
        M_old = M;  M = M/betarho;
        % else
        %    rho = rho*betarho;
    end
    
    %A=P*F*P+rho*Q; % A=F+rho*Q;
    b_old=b;     b=b0-G'*mi;
    nout = nout + 1;
    if (ncg>=ncg_max) break; end;
end
geom.plot(u(1:n));

[nout ncg ne np]
figure; semilogy(Vgp, 'r'); hold on; semilogy(Vnarus, 'g'); semilogy(abs(VLagIn), 'b'); semilogy(VCrit, 'k');
figure; semilogy(abs(VLag), 'r'); hold on; semilogy(VM, 'b');