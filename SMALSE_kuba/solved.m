function [flg] = solved(nrmg,x,ctx)
% SMALXE convergence test
  nrmgeq = norm(ctx.B*x-ctx.c);
  innertol = min(ctx.M*nrmgeq,ctx.eta);
  %innertol = min(ctx.M*ctx.nrmgeq,ctx.eta);
  %innertol = max(ctx.M*ctx.nrmgeq,innertol);
  %innertol = ctx.M*ctx.nrmgeq;
  flg = 0; % did not converge
  if nrmg <= ctx.outertol && nrmgeq <= ctx.outertol
    flg = 2; % happy breakdown
  else if nrmg <= innertol
    flg = 1; % innersolver convergence
  end
  %printf("%e %e %e %e\n",nrmg,nrmgeq,innertol,ctx.outertol);
end

