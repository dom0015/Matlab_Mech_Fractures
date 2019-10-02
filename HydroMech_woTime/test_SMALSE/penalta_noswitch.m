function [u,i,d_vec,un_vec,geom_vec,rho] = penalta_noswitch(A,b,B_func,d_min,max_it,d_eps,rho_step,rho_step_limits,stagnate_count,alpha)
rho_max=1e0;
rho_min=1e-20;
B=B_func(0*b);
rho=sparse(ones(size(B,1),1));

c=-d_min*ones(size(B,1),1);
H_P=2*B'*diag(rho)*B;
u=0*b;
u_old=u;
grad_P=H_P*u-(2*c'*diag(rho)*B)';
d_vec=zeros(max_it,1);
un_vec=zeros(max_it,1);
geom_vec=zeros(max_it,1);

idx_on_old=0*c;
idx_on_changes=0*idx_on_old;
B_old=B;
l_vypis=0;
for i=1:max_it
    
    u_n=(A+H_P)\(A*u+grad_P-b);   
    u=u-alpha*u_n;
    
    d=B*u-c;
    d(d<0)=0;
    norm_d=max(d);
    
    B=B_func(u);
    H_P=2*B'*diag(rho)*B;
    grad_P=H_P*u-(2*c'*diag(rho)*B)';
   
    d=B*u-c;
    idx_on=d>0;
    
    res_geom=full(max(max(abs(B_old-B))));
    norm_un=norm(u_n);
    d_vec(i)=norm_d;
    un_vec(i)=norm_un;
    geom_vec(i)=res_geom;
    
    if d_vec(i)<d_eps  && norm_un<d_eps && res_geom<d_eps
        fprintf('Converged: penalty=%d in it=%d, norm(d)=%d\n',max(full(rho)),i,d_vec(i));
        break
    end
    
    if  max(idx_on_changes)>stagnate_count %d_vec(i)>stagnate_ratio*min(d_vec(1:i)) ||
        u=u_old;
        B=B_func(u);
        d=B*u-c;
        idx_on=d>eps;
        fprintf('Penalty stagnated: penalty=%.1f/%.1f in it=%d, norm(d)=%d, stagnated norm(d)=%d\n',full(log10(max(rho))),full(log10(min(rho))),i,max(d(idx_on)),d_vec(i));
        break;
    end
    
    l_vypis = my_print(sprintf('It. %d: Penalty= %.1f/%.1f, norm(d)=%d, geom=%d, un=%d, itch=%d \n',i,full(log10(max(rho))),full(log10(min(rho))),max(d(idx_on)),res_geom,norm(u_n),max(idx_on_changes)),l_vypis );
    
    idx_on_changes=idx_on_changes+~(idx_on_old==idx_on);

    
    rho(idx_on)=rho(idx_on)*rho_step;
    rho(~idx_on)=rho(~idx_on)/rho_step;
    rho=min(rho,rho_max);
    rho=max(rho,rho_min);
    
    if norm_un<10^(-log10(rho_max)) && res_geom<max(10^(-log10(rho_max)))
        rho_max=rho_step_limits*rho_max;
        idx_on_changes(:)=0;
        rho_min=rho_min/rho_step_limits;
    end
    idx_on_old=idx_on;
    u_old=u;
    B_old=B;
end
d_vec=d_vec(1:i);
un_vec=un_vec(1:i);
geom_vec=geom_vec(1:i);
if i==max_it
    fprintf('Maximum iterations reached: penalty=%d in it=%d, norm(d)=%d\n',full(log10(max(rho))),i,d_vec(end));
end
rho=full(rho);
end

