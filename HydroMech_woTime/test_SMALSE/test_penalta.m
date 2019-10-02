%%
n_node=101; % 5*k+1;
[A,b,B,geom] = SMALSE_precomp(n_node);

d_eps=1e-8;
d_min=0;
max_it=1000;
rho_step=3;
rho_step_limits=10;
stagnate_count=10;
alpha=1;

[u,niter,d_vec,un_vec,geom_vec,rho] = penalta_noswitch(A,b,geom.recalculate_B,d_min,max_it,d_eps,rho_step,rho_step_limits,stagnate_count,alpha);

figure
plot(d_vec,'r');
set(gca,'YScale','log')
hold on
plot(un_vec,'b');
plot(geom_vec,'g');

figure
plot(rho)
set(gca,'YScale','log')

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