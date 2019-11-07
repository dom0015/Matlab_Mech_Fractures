function [flg,nrmg,nrmgeq] = solved_fnc1(nrmg,x,ctx)
% SMALXE convergence test
nrmgeq = norm(ctx.B(x)-ctx.c);
[flg,nrmg,nrmgeq]=solved_fnc(nrmg,nrmgeq,ctx);
end
