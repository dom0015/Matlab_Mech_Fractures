%%
n_node=101; % 5*k+1;
[A,b,B,geom] = SMALSE_precomp(n_node);

d_eps=1e-12;
d_min=0;
max_it=10000;
rho_step=10;


[u_0] = penalta(A,b,geom.recalculate_B,d_min,max_it,1e-8,2);

[u,niter] = bariera(A,b,geom.recalculate_B,d_min,max_it,rho_step,d_eps,geom.plot,u_0);

geom.plot(u);
% 
% 
% 
% 
% 
% 
% 
% 
% M=[A B(idx_on,:)';B(idx_on,:) zeros(sum(idx_on))];
% b_lambda=[b;c(idx_on)];
% u_idx=1:length(u);
% lambda_idx=(length(u)+1):length(b_lambda);
% u_lambda=M\b_lambda;
% geom.plot(u_lambda(u_idx));
% u_lambda(lambda_idx)'
% d=B*u_lambda(u_idx);
% (d(d>c)-c((d>c)))'