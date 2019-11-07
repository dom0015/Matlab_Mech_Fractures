function [x,flg,k,iter,errors] = mprgp_fnc(A,x,b,ctx,maxit,mprgpctx)
% solve .5x'Ax -x'b s.t. x >= lb
%
% ctx = covergence ctx
% % ctx.outertol = rtol*norm(b0);
% % ctx.B = constraint matrix
% % ctx.c = constraint vec
% % ctx.M = M;
% % ctx.eta = eta;
% % ctx.rhochanged = need to recompute abar
% % ctx.nrmgeq = ||Bx-c||
%
% mprgpctx
% % mprgpctx.lb = bound
% % mprgpctx.abarmult = (0,2], length of expansion step = abarmult*||A^{-1}||
% % mprgpctx.propConst = proportioning const.
% % mprgpctx.settol = 10*eps, splitting tolerance
% % mprgpctx.infeastol = 0, infeasibility tolerance

k = [0 0 0 0]; %Hessian, CG, Exp, Prop
iter = 0;
if  ~exist("mrpgpctx.abar","var") %||ctx.rhochanged
    if ctx.rhochanged
        ctx.rhochanged = false;
    end
    mprgpctx.abar = mprgpctx.abarmult/ctx.rho;
    %mprgpctx.abar = mprgpctx.abarmult/mprgpctx.maxeig;
end
x = max(x,mprgpctx.lb); % project to feasible set
g = A(x) - b; k = k + [1 0 0 0];
[gf,gc,gr] = mprgpSplit(x,g,mprgpctx);
p = gf;
[flg,nrmg,nrmgeq] = solved_fnc1(norm(gf+gc),x,ctx);
errors=zeros(maxit,2);
errors(1,:)=[nrmg nrmgeq];
while ~flg && iter < maxit
    if gc'*gc <= mprgpctx.propConst^2 *gr'*gf % proportional
        Ap = A(p); k = k + [1 0 0 0];
        pAp = p'*Ap;
        acg = g'*p/pAp;
        bounddiff = x - mprgpctx.lb;
%       afeas = Inf;
%         for i=1:length(x)
%             %if (p(i) > mprgpctx.settol)
%             if (p(i) > 0)
%                 af = bounddiff(i)/p(i);
%                 if af < afeas
%                     afeas = af;
%                 end
%             end
%         end
        af=bounddiff./p;
        af(p <= 0)=inf;
        afeas=min(af);
        
        
        
        if acg <= afeas % cg % step = 'c';
            
            x = x - acg*p;
            g = g - acg*Ap;
            [gf,gc,gr] = mprgpSplit(x,g,mprgpctx);
            bcg = gf'*Ap/pAp;
            p = gf - bcg*p;
            k = k + [0 1 0 0];
        else % expansion %step = 'e';
            x = x - afeas*p;
            g = g - afeas*Ap;
            [gf,gc,gr] = mprgpSplit(x,g,mprgpctx);
            %x = x - mprgpctx.abar*gr;
            x = x - mprgpctx.abar*gf;
            x = max(x,mprgpctx.lb);
            g = A(x) - b; k = k + [1 0 1 0];
            [gf,gc,gr] = mprgpSplit(x,g,mprgpctx);
            p = gf;
        end
    else % proportioning step = 'p';
        p = gc;
        Ap = A(p); k = k + [1 0 0 1];
        pAp = p'*Ap;
        acg = g'*p/pAp;
        x = x - acg*p;
        g = g - acg*Ap;
        [gf,gc,gr] = mprgpSplit(x,g,mprgpctx);
        p = gf;
    end
    %flg = solved_fnc(norm(gf+gc),x,ctx);
    [flg,nrmg,nrmgeq] = solved_fnc(norm(gf+gc),nrmgeq,ctx);
    iter = iter + 1;
    errors(iter,:)=[nrmg nrmgeq];
end
errors=errors(1:max(iter,1),:);
if flg==0
    flg=1;
end
end


