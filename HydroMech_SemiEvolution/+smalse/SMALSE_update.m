function [u,mi,update_data,nout,ncg] = SMALSE_update(F,b0,G,c,update_G,u0,mi0,update_data,idx_no_bounds,rel,tol_to_update,rho0,betarho,Gama,M_start,maxiter_cg,type,printres)
%% SMALSE  with updates of geometry for dual formulation
% includes initial guess, adaptive change of rho/M
% INPUTS -----------------------------------------------
% F - matrix from minimization u'Fu-u'b0
% b0 - rhs from minimization
% G - equality constrain Gu=0
% c - lower bounds u>=c >>>>> u=lambda_ker,I, c=-lambda_Im,I
% u0 - original guess of the solution
% update_Geom - [F,b0,G,c,update_data,geom_res] = update_G(u,update_data);
%          - fnc for updating geometry
% update_data - changing data on update struct(lambda_ImGt, B_i)
% idx_no_bounds - indexes of nodes without lower bound (c=-inf)
% epsr - precision for norm(gp) and norm(Gu)
% epsrGeom - precision for geometry change (geom_res of update_Geom)
% epsupd - maximum value of abs(Lag-LagOld) for allowed update of params
% tol_to_update - minimum precision of norm(Gu) for geometry update
% rho0 - starting value of rho (lFl=lG'Gl=1)
% betarho - increment of rho=betarho*rho / M=M/betarho
% Gama - gama parameter (gc'*gc <= Gama^2*gr'*gf)
% M_start - starting value of M (M=M_start*lAl)
% maxiter_cg - maximum of cg steps
% type - type of update ('M'-change only M, 'rho' - change only rho, 
%        'm' - change M and rho depending on gp/Gu decrease)
% printres - true/false - printig intermediate results and graphs
%% -------------------------------------------------------

c(idx_no_bounds) = -Inf*ones(length(idx_no_bounds),1);
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

if ~isempty(update_G)
    geom_res=-1;
else
    geom_res=0;
end
Vgeom_res=[];

b=b0-G'*mi;


A=F+rho*Q;
lAl=eigs(A, 1);
alfa=1/lAl;
M=M_start*lAl;

Lag=0.5*u'*A*u-b'*u;
Lag_old=Lag;
gp=ones(length(u),1);
norm_b0=norm(b0);

l_vypis=0;

Gu_old=Gu;
gp_old=gp;
updatednow=false;
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
        if norm(gp)<crit || (norm(gp)<epsr && norm(Gu)<epsr && abs(geom_res)<rel)
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
    if Lag - Lag_old < Increment && (abs(Lag - Lag_old)>epsr*norm_b0 ) && ~updatednow
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
            if norm(gp)/norm(gp_old)<0.75*norm(Gu)/norm(Gu_old)
                rho=rho*betarho;
                A=F+rho*Q;
                lAl=max(1,rho);
                alfa=1/lAl;
            else
                M = M/betarho;
            end
        end
    end
    
    if (norm(gp)<epsr && norm(Gu)<epsr)% && abs(geom_res)<rel
        if printres
        fprintf('Converged: norm(gp)=%d  norm(Gu)=%d.\n',norm(gp),norm(Gu));
        end
        break
    end
    updatednow=false;
    if ~isempty(update_G)
        if norm(Gu)<tol_to_update
            tol_to_update=max(norm(Gu)/10,epsr);
            [F,b0,G,c,update_data,geom_res,~,alpha] = update_G(u,update_data);
            F=F/lFl_*lQl_;
            b0=b0/lFl_*lQl_;
            
            u=max(u,c);
            Gu=G*u;
            mi=-alpha/lFl_*lQl_;
            mi=mi+rho*Gu;
            J = (u > c);
            Q=G'*G;
            A=F+rho*Q;
            updatednow=true;
            Vgeom_res=[Vgeom_res geom_res];
        else
            Vgeom_res=[Vgeom_res -1];
        end
    else
        Vgeom_res=[Vgeom_res -1];
    end
    
    
    Gu_old=Gu;
    gp_old=gp;
    
    b=b0-G'*mi;
    nout = nout + 1;
    if ~isempty(Vnarus)
    if printres
        l_vypis = smalse.my_print( sprintf('Outer it=%d CG iter=%d norm(gp)=%d norm(G*u)=%d, M=%d, rho=%d\n',nout, ncg,Vgp(end),Vnarus(end),M,rho),l_vypis );
    end
    end
    if (ncg>=maxiter_cg)
        fprintf('Max. iter. ');
        fprintf('Outer it=%d CG iter=%d norm(gp)=%d norm(G*u)=%d, M=%d, rho=%d\n',nout, ncg,Vgp(end),Vnarus(end),M,rho)
        break
    end
end
Vgp=[Vgp norm(gp)];
Vnarus=[Vnarus norm(Gu)];
fprintf('Outer it=%d CG iter=%d norm(gp)=%d norm(G*u)=%d, M=%d, rho=%d\n',nout, ncg,Vgp(end),Vnarus(end),M,rho)

if printres
fprintf('Number of iterations: nout=%d ncg=%d ne=%d np=%d\n', nout, ncg, ne, np);
end
if printres
    figure(101);
    semilogy(Vgp, 'r'); hold on;
    semilogy(Vnarus, 'g');
    semilogy(abs(VLagIn), 'b');
    semilogy(VCrit, 'k');
    semilogy(Vcg,abs(VLag), 'r+:');
    semilogy(Vcg,VM, 'b+-');
    semilogy(Vcg,Vrho, 'k+-');
    semilogy(Vcg,VIncrement, 'g+-');
    
    if (ncg>=maxiter_cg)
%         semilogy(Vcg,Vgeom_res, 'k*:');
    else
        semilogy(Vcg(1:end-1),Vgeom_res, 'k*:');
    end
    semilogy(Vexpng,Vexpgp, 'b+','MarkerSize',10);
    title('Convergence in CG iterations');
    xlabel('CG iter');
    ylabel('value');
    legend({'gp','narus','abs(LagIn)','Crit','abs(Lag)','M','rho','increment','Bi-Biold','expanze'});
    hold off
end

% global S_iter
% S_iter=[S_iter ncg];

% is_contact = false;
% temp = u(end-38:end)-c(end-38:end);
% if max(temp)>0
%     is_contact = true;
% end
% figure(300); plot(temp); hold on; title('multiplicators')
end

