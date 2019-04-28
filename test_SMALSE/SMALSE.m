function [u] = SMALSE(F,b0,G,c,idx_no_bounds,rel,rho0,betarho,Gama,M_start,maxiter_cg)
% SMALSE  David
% MPRGP
% INPUTS -----------------------------------------------
% F -
% b0 -
% G -
%-------------------------------------------------------

lFl=rayleigh_est(F, 1000,1e-2);

epsr = rel*norm(b0);

ncg = 0; 
ne = 0; 
np = 0; 
nout = 0;

VLag=[];
VLagIn=[];
Vnarus=[];
Vgp=[];
VM=[];
VCrit=[];

c(idx_no_bounds) = -Inf*ones(length(idx_no_bounds),1);


    u0=zeros(length(c),1);

u=max(u0,c+eps);


J = (u > c);
rho = rho0* lFl;
Gu=G*u;
mi=zeros(size(G,1),1);  
mi=mi+rho*Gu;
b=b0-G'*mi;
Q=G'*G;
A=F+rho*Q; % A=P*F*P+rho*Q;
lAl=rayleigh_est(A, 1000,1e-2);
alfa=1/lAl;
M=M_start*lAl;

Lag=0.5*u'*A*u-b'*u;
gp=ones(length(u),1);
norm_b0=norm(b0);

l_vypis=0;

while (1)
    if (norm(gp)<epsr && norm(Gu)<epsr)
        fprintf('Converged: norm(gp)=%d  norm(Gu)=%d.\n',norm(gp),norm(Gu));
        break
    end
    
    % Gradient splitting of g=-r, gf=gradient free(fi), gr=gradient free
    % reduced (fired), gc=gradient chopped or cut (beta), gp=gf+gc=projected grad.
    % previously used gr = min(lAl*J.*(u-c), gf);
    
%     if nargin>=11
%         if ~isempty(update_G)
%             G=update_G(u); 
%             Q=G'*G;
%             A=F+rho*Q;
%         end
%     end  
    
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
    VLag=[VLag Lag];
    VM=[VM M];
    Lag_old=Lag;
    Lag=0.5*u'*A*u-b'*u;
    
    Increment=0.5*rho*crit^2;
    
    if Lag - Lag_old < Increment

            M = M/betarho;

    end
    
    b=b0-G'*mi;
    nout = nout + 1;
    
    l_vypis = my_print( sprintf('Outer it=%d CG iter=%d norm(gp)=%d norm(G*u)=%d\n',nout, ncg,Vgp(end),Vnarus(end)),l_vypis );
    
    if (ncg>=maxiter_cg)
        fprintf('Maximum iterations reached.\n');
        break
    end
end
fprintf('Number of iterations: nout=%d ncg=%d ne=%d np=%d\n', nout, ncg, ne, np);
figure; 
semilogy(Vgp, 'r'); hold on;  
semilogy(Vnarus, 'g'); 
semilogy(abs(VLagIn), 'b'); 
semilogy(VCrit, 'k');
title('Convergence in CG iterations')
xlabel('CG iter')
ylabel('value')
legend({'gp','narus','abs(LagIn)','Crit'})

figure; 
semilogy(abs(VLag), 'r'); hold on; 
semilogy(VM, 'b');
title('Convergence in outer iter')
xlabel('Outer iter')
ylabel('value')
legend({'abs(Lag)','M'})
end

