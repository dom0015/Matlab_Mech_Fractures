function [u,update_data,nout,ncg] = SMALSE_update(F,b0,G,c,update_G,update_data,idx_no_bounds,rel,tol_to_update,rho0,betarho,Gama,M_start,maxiter_cg,type,printres)
% SMALSE  David
% MPRGP
% INPUTS -----------------------------------------------
% F -
% b0 -
% G -
%-------------------------------------------------------


c(idx_no_bounds) = -Inf*ones(length(idx_no_bounds),1);


u0=zeros(length(c),1);

u=max(u0,c);

% if ~isempty(update_G)
%     
%     [F,b0,G,c,update_data] = update_G(u,update_data);
%     
% end


Q=G'*G;
lFl_=eigs(F, 1);
lQl_=eigs(Q, 1);

F=F/lFl_*lQl_;
b0=b0/lFl_*lQl_;
lFl=eigs(F, 1);
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

J = (u > c);
rho = rho0* lFl;
Gu=G*u;
mi=zeros(size(G,1),1);
mi=mi+rho*Gu;
b=b0-G'*mi;


A=F+rho*Q;
lAl=eigs(A, 1);
alfa=2/lAl;
M=M_start*lAl;

Lag=0.5*u'*A*u-b'*u;
Lag_old=Lag;
gp=ones(length(u),1);
norm_b0=norm(b0);

l_vypis=0;

Gu_old=Gu;
gp_old=gp;
tag=false;
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
    p_all=[];
    g_all=[];
    while (1)
        if ncg==1362
        end
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
%                 if ~isempty(g_all)
%                     gp=gp-g_all(:,max(1,end-200):end)*(g_all(:,max(1,end-200):end)'*gp);
% %                     gf = J.*g;
% %                     gc = min((~J).*g,0);
% %                     gp = gf + gc;
%                 end
%                 g_all=[g_all gp/norm(gp)];
                if size(g_all,2)>10
                end
                beta=gf'*Ap/pAp;
                p=gf-beta*p;
%                 if ~isempty(p_all)
%                     p=p-p_all(:,max(1,end-200):end)*(p_all(:,max(1,end-200):end)'*(A*p));
%                     p=p-p_all(:,max(1,end-200):end)*(p_all(:,max(1,end-200):end)'*(A*p));
%                     p=p-p_all(:,max(1,end-200):end)*(p_all(:,max(1,end-200):end)'*(A*p));
%                 end
%                 p_all=[p_all p/sqrt((p'*A*p))];
            else
                p_all=[];
                g_all=[];
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
    
    if Lag - Lag_old < Increment && (abs(Lag - Lag_old)>epsr )
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
            if norm(gp)/norm(gp_old)<min(M,1)*norm(Gu)/norm(Gu_old)
                rho=rho*betarho;
                A=F+rho*Q;
                lAl=max(1,rho);
                alfa=1/lAl;
            else
                M = M/betarho;
            end
        end
    end
    
    if (norm(gp)<epsr && norm(Gu)<epsr) && abs(geom_res)<rel
        fprintf('Converged: norm(gp)=%d  norm(Gu)=%d.\n',norm(gp),norm(Gu));
        break
    end
    
    if ~isempty(update_G)
        geom_res=-1;
        if norm(Gu)<tol_to_update %&& ~tag
            tol_to_update=max(norm(Gu)/10,epsr);
      %      betarho=1;
            tag=true;
            
            [F,b0,G,c,update_data,geom_res] = update_G(u,update_data);
            F=F/lFl_*lQl_;
            b0=b0/lFl_*lQl_;
            
            u=max(u,c);
            Gu=G*u;
%             J = (u > c+norm(Gu)^(0.8));
%             u(~J)=c(~J);
%             Gu=G*u;
%             tres=(b0-F*u);
%             mi=(G(:,J)')\(tres(J));
%             
            J = (u > c);
            Q=G'*G;
            A=F+rho*Q;
        end
    end
    
    Vgeom_res=[Vgeom_res geom_res];
    Gu_old=Gu;
    gp_old=gp;
    
    b=b0-G'*mi;
    nout = nout + 1;
    if printres
        l_vypis = my_print( sprintf('Outer it=%d CG iter=%d norm(gp)=%d norm(G*u)=%d, M=%d, rho=%d\n',nout, ncg,Vgp(end),Vnarus(end),M,rho),l_vypis );
    end
    if (ncg>=maxiter_cg)
        fprintf('Maximum iterations reached.\n');
        break
    end
end
fprintf('Number of iterations: nout=%d ncg=%d ne=%d np=%d\n', nout, ncg, ne, np);
if printres
    figure;
    semilogy(Vgp, 'r'); hold on;
    semilogy(Vnarus, 'g');
    semilogy(abs(VLagIn), 'b');
    semilogy(VCrit, 'k');
    semilogy(Vcg,abs(VLag), 'r+:');
    semilogy(Vcg,VM, 'b+-');
    semilogy(Vcg,Vrho, 'k+-');
    semilogy(Vcg,VIncrement, 'g+-');
    
    if (ncg>=maxiter_cg)
        semilogy(Vcg,Vgeom_res, 'k*:');
    else
        semilogy(Vcg(1:end-1),Vgeom_res, 'k*:');
    end
    semilogy(Vexpng,Vexpgp, 'b+','MarkerSize',10);
    title('Convergence in CG iterations');
    xlabel('CG iter');
    ylabel('value');
    legend({'gp','narus','abs(LagIn)','Crit','abs(Lag)','M','rho','increment','Bi-Biold','expanze'});
end
end

