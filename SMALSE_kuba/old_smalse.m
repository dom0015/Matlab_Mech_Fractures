function [u,mi,nout,ncg] = old_smalse(F,b0,G,c,u0,mi0,rel,rho0,betarho,Gama,M_start,maxiter_cg,type,printres)
%% SMALSE  with updates of geometry for dual formulation
% includes initial guess, adaptive change of rho/M
% INPUTS -----------------------------------------------
% F - matrix from minimization u'Fu-u'b0
% b0 - rhs from minimization
% G - equality constrain Gu=0
% c - lower bounds u>=c >>>>> u=lambda_ker,I, c=-lambda_Im,I
% u0 - original guess of the solution
% rel - precision for norm(gp) and norm(Gu)
% rho0 - starting value of rho (lFl=lG'Gl=1)
% betarho - increment of rho=betarho*rho / M=M/betarho
% Gama - gama parameter (gc'*gc <= Gama^2*gr'*gf)
% M_start - starting value of M (M=M_start*lAl)
% maxiter_cg - maximum of cg steps
% type - type of update ('M'-change only M, 'rho' - change only rho,
%        'm' - change M and rho depending on gp/Gu decrease)
% printres - true/false - printig intermediate results and graphs
%% -------------------------------------------------------

Q=G'*G;
lFl_=eigs(F, 1);
lQl_=1;

F=F/lFl_*lQl_;
b0=b0/lFl_*lQl_;
lFl=eigs(F, 1);
rho = rho0* lFl;

if isempty(u0)
    u0=zeros(length(c),1);
    u=max(u0,c);
    J = (u > c);
    Gu=G*u;
    mi=zeros(size(G,1),1);
    mi=mi+rho*Gu;
else
    u=max(u0,c);
    Gu=G*u;
    mi=mi0;
    mi=mi+rho*Gu;
    J = (u > c);
end

epsr = rel*norm(b0);

ncg = 0;
ne = 0;
np = 0;
nout = 0;

VLag=[];
Vcg=[];
VLagIn=[];
Vnarus=[];
Vgp=[];
VM=[];
Vrho=[];
VCrit=[];
VIncrement=[];
Vexpng=[];
Vexpgp=[];

b=b0-G'*mi;

A=F+rho*Q;
lAl=eigs(A, 1);
alfa=1/lAl;
M=M_start*lAl;

Lag=0.5*u'*A*u-b'*u;
Lag_old=Lag;
gp=ones(length(u),1);
norm_b0=norm(b0);

Gu_old=Gu;
gp_old=gp;

while (1)
    % Gradient splitting of g=-r, gf=gradient free(fi), gr=gradient free
    % reduced (fired), gc=gradient chopped or cut (beta), gp=gf+gc=projected grad.
    % previously used gr = min(lAl*J.*(u-c), gf);
    g = A*u - b;
    gf = J.*g;
    gc = min((~J).*g,0);
    gr = min(J.*(u-c)/alfa, gf);
    gp = gf + gc;
    p=gf;
    Gu=G*u;
    while (1)
        crit=min(M*norm(Gu),norm_b0);
        if norm(gp)<crit || (norm(gp)<epsr && norm(Gu)<epsr)
            break
        end
        VCrit=[VCrit crit];
        
        if gc'*gc <= Gama^2*gr'*gf
            % Proportional iteration. Trial conjugate gradient step.
            Ap=A*p;
            rtp=g'*p;
            pAp=p'*Ap;
            acg=rtp/pAp;
            yy=u-acg*p;
            ncg=ncg+1;
            if all(yy>=c)
                %Conjugate gradient step
                u=yy;
                g=g-acg*Ap;
                
                gf = J.*g;
                gc = min((~J).*g,0);
                gr = min(lAl*J.*(u-c), gf);
                gp = gf + gc;
                beta=gf'*Ap/pAp;
                p=gf-beta*p;
            else
                %partial conjugate gradient step
                if max((J.*p)./(J.*(u-c)+(~J)))~=0
                    a=1/max((J.*p)./(J.*(u-c)+(~J)));
                else
                    a=0;
                end
                u=max(u-a*p,c);
                J=(u>c);
                g=g-a*Ap;
                gf = J.*g;
                gr = min(lAl*J.*(u-c), gf);
                %expansion step
                u=max(u-alfa*gr,c);
                J=(u>c);
                g=A*u-b;
                gf = J.*g;
                gc = min((~J).*g,0);
                gr = min(lAl*J.*(u-c), gf);
                gp = gf + gc;
                p=gf;
                ne=ne+1;
                Vexpng=[Vexpng ncg];
                Vexpgp=[Vexpgp norm(gp)];
            end
        else
            %Proportioning step
            Ap=A*gc;
            acg=(gc'*g)/(gc'*Ap);
            u=u-acg*gc;
            J=(u>c);
            g=g-acg*Ap;
            gf = J.*g;
            gc = min((~J).*g,0);
            gr = min(lAl*J.*(u-c), gf);
            gp = gf + gc;
            p=gf;
            np=np+1;
        end
        
        Gu=G*u;
        VLagIn=[VLagIn 0.5*u'*(g-b)];
        Vnarus=[Vnarus norm(Gu)];
        Vgp=[Vgp norm(gp)];
        if (ncg>=maxiter_cg)
            break
        end
    end
    
    mi=mi+rho*Gu;
    VLag=[VLag abs(Lag - Lag_old)];
    VM=[VM M];
    Vrho=[Vrho rho];
    Vcg=[Vcg ncg];
    Lag_old=Lag;
    Lag=(0.5*u'*A-b')*u;
    
    Increment=0.5*rho*norm(Gu)^2;
    VIncrement=[VIncrement Increment];
    
    % update M provided the increase of the Lagrangian in not sufficient
    if Lag - Lag_old < Increment && (abs(Lag - Lag_old)>epsr*norm_b0 )
        if type=='M'
            M = M/betarho;
        end
        if type=='r'
            rho=rho*betarho;
            A=F+rho*Q;
            lAl=max(1,rho);
            alfa=1/lAl;
        end
        if type=='m'
            if norm(gp)/norm(gp_old)<norm(Gu)/norm(Gu_old)
                rho=rho*betarho;
                A=F+rho*Q;
                lAl=max(1,rho);
                alfa=1/lAl;
            else
                M = M/betarho;
            end
        end
    end
    
    if (norm(gp)<epsr && norm(Gu)<epsr)
        if printres
            fprintf('Converged: norm(gp)=%d  norm(Gu)=%d.\n',norm(gp),norm(Gu));
        end
        break
    end
    
    Gu_old=Gu;
    gp_old=gp;
    
    b=b0-G'*mi;
    nout = nout + 1;
    if ~isempty(Vnarus)
        if printres
            fprintf('Outer it=%d CG iter=%d norm(gp)=%d norm(G*u)=%d, M=%d, rho=%d\n',nout, ncg,Vgp(end),Vnarus(end),M,rho);
        end
    end
    if (ncg>=maxiter_cg)
        fprintf('Max. iter. ');
        fprintf('Outer it=%d CG iter=%d norm(gp)=%d norm(G*u)=%d, M=%d, rho=%d\n',nout, ncg,Vgp(end),Vnarus(end),M,rho)
        break
    end
end
fprintf('Outer it=%d CG iter=%d norm(gp)=%d norm(G*u)=%d, M=%d, rho=%d\n',nout, ncg,Vgp(end),Vnarus(end),M,rho)

if printres
    fprintf('Number of iterations: nout=%d ncg=%d ne=%d np=%d\n', nout, ncg, ne, np);
end

end

