function [x,flg,k,mu] = smalxe(A,b,B,c,x0,rho,rho_update,M,M_update,eta,rtol,maxit,innersolver,varargin)
% Solves QP s.t Bx = c, where QP is solved by inner solver
% smalxe(A,b,B,c,x0,10,1,10,10,1e-1,1e-6,innersolver,...)
% innersolver = @cgCommon
% ... innersolver arguments

% initialize
x = x0;
b0 = b;
maxeig = eigs(A,1);
rho = rho*maxeig;
M = M*maxeig;
eta = eta*norm(b0);
%eta = 1e12;
Apenalized = A + rho*B'*B;
mu = zeros(size(B,1),1); 
%mu=mu+rho*(B*x-c);
b = b0 - B'*(mu + rho*c);
lag = 0.5*x'*Apenalized*x - b'*x;

% convergence ctx
ctx.outertol = rtol*norm(b0);
ctx.B = B;
ctx.c = c;
ctx.M = M;
ctx.eta = eta;
ctx.rhochanged = true;
ctx.nrmgeq = norm(b0);

kinner = 0;
kouter = 0;
while kouter < maxit
  % minimize
  [x,flg,k] = innersolver(Apenalized,x,b,ctx,varargin{:});
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
  lag = 0.5*x'*Apenalized*x -b'*x;
  increment = 0.5*rho*norm(B*x - c)^2;
  
  % update multiplicator
  mu = mu + rho*(B*x - c);

  % update M, rho
  if lag - lag_old < increment   
    M = M/M_update;
    ctx.M = M;
    % rho variant
%   rho = rho*rho_update; 
%   Apenalized = A + rho*B'*B;
%   ctx->rhochanged = 1;
  end

  % update inner rhs
  b = b0 - B'*(mu + rho*c); 
  fprintf("norm Gx: %d\n",norm(B*x - c));
%ctx.nrmgeq = norm(b);
end

flg = 0;

