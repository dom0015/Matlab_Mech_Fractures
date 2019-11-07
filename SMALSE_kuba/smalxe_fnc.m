function [x,flg,k,mu,all_errors,all_iters] = smalxe_fnc(A,b,B,Bt,c,x0,mu0,rho,rho_update,M,M_update,eta,rtol,maxit,innersolver,type,varargin)
% Solves QP s.t Bx = c, where QP is solved by inner solver
% smalxe(A,b,B,c,x0,10,1,10,10,1e-1,1e-6,innersolver,...)
% innersolver = @cgCommon
% ... innersolver arguments

% initialize
x = x0;
b0 = b;
maxeig = eigs(A,length(b),1);
rho = rho*maxeig;
M = M*maxeig;
eta = eta*norm(b0);
%eta = 1e12;
Apenalized =@(x) A(x) + rho*Bt(B(x));
mu = mu0;

% if isempty(x0)
%     u0=zeros(length(c),1);
%     u=max(u0,c);
%     J = (u > c);
%     Gu=G*u;
%     mi=zeros(size(G,1),1);
%     mi=mi+rho*Gu;
% else
%     u=max(u0,c);
%     Gu=G*u;
%     mi=mi0;
%     mi=mi+rho*Gu;
%     J = (u > c);
% end

%mu=mu+rho*(B*x-c);
b = b0 - Bt(mu + rho*c);
lag = 0.5*x'*Apenalized(x) - b'*x;

% convergence ctx
ctx.outertol = rtol*norm(b0);
ctx.B = B;
ctx.Bt = Bt;
ctx.c = c;
ctx.M = M;
ctx.eta = eta;
ctx.rhochanged = true;
ctx.nrmgeq = norm(b0);
ctx.rho=rho;

kinner = 0;
kouter = 0;
all_errors=[];
all_iters=[];
while kouter < maxit
    % minimize
    [x,flg,k,~,error] = innersolver(Apenalized,x,b,ctx,varargin{:});
    all_errors=[all_errors;error];
    all_iters=[all_iters;k];
    current_notmGU=norm(B(x) - c);
    fprintf("norm Gp: %d Gx: %d\n",all_errors(end,1),current_notmGU);
    fprintf("M = %d, rho = %d\n",M,rho)
    if ~kinner
        kinner = k;
    else
        kinner = kinner + k;
    end
    kouter = kouter + 1;
    k = [kouter, kinner];
    
    % outer stopping criteria
    if ~flg
        fprintf("Inner solver diverged\n");
        return
    elseif flg == 2
        flg = 1;
        return
    end
    
    %% compute Lagrangian
    lag_old = lag;
    lag = 0.5*x'*Apenalized(x) -b'*x;
    increment = 0.5*rho*norm(B(x) - c)^2;
    
    % update multiplicator
    mu = mu + rho*(B(x) - c);
    
    % update M, rho
    if lag - lag_old < increment
        if type=='M'
            M = M/M_update;
            ctx.M = M;
        end
        if type=='r'
            rho=rho*rho_update;
            Apenalized =@(x) A(x) + rho*Bt(B(x));
            ctx.rho=rho;
        end
        if type=='m'
            if error(end,1)/error(1,1)<current_notmGU/error(1,2)
                rho=rho*rho_update;
                Apenalized =@(x) A(x) + rho*Bt(B(x));
                ctx.rho=rho;
            else
                M = M/M_update;
                ctx.M = M;
            end
        end
    end
    
    
    
    % update inner rhs
    b = b0 - Bt(mu + rho*c);
    %ctx.nrmgeq = norm(b);
end

flg = 0;

