function [x,flg,k,iter] = old_mprgp(A,x,b,ctx,maxit,mprgpctx)

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

end
