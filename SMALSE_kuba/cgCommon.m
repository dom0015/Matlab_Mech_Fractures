function [x,flg,k] = cgCommon(A,x,b,ctx,maxit)
    k = 0;
    r = b - A*x;
    p = r;
    alpha_n = r'*r;
    flg = solved(sqrt(alpha_n),x,ctx);
    while ~flg && k < maxit
        s = A*p;
        alpha = alpha_n/(p'*s);
        x = x + alpha*p;
        r = r - alpha*s;
        beta_d = alpha_n;
        alpha_n = r'*r;
        beta = alpha_n/beta_d;
        p = r + beta*p;
        k = k+1;
        flg = solved(sqrt(alpha_n),x,ctx);
    end
end
