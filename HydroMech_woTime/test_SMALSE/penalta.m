function [u,i,d_vec] = penalta(A,b,B_func,d_min,max_it,d_eps,rho_step)
stagnate_ratio=10;

B=B_func(0*b);
rho=1/size(B,1);

H_P=sparse(size(A,1),size(A,2));
u=0*b;
grad_P=H_P*u;
d_vec=zeros(max_it,1);
c=-d_min*ones(size(B,1),1);
B=0*B;
B_old=B;
l_vypis=0;
for i=1:max_it
    u=u-(A+rho*H_P)\(A*u+rho*grad_P-b);
    B=B_func(u);
    d=B*u-c;
    idx_on=d>eps;
    H_P=2*B(idx_on,:)'*B(idx_on,:);
    grad_P=H_P*u-(2*c(idx_on)'*B(idx_on,:))';
    norm_d=max(d(idx_on));
    d_vec(i)=norm_d;
    
    res_geom=full(max(max(abs(B_old-B))));
    
    if d_vec(i)<d_eps
        fprintf('Converged: penalty=%d in it=%d, norm(d)=%d\n',rho,i,d_vec(i));
        break
    end
    
    if d_vec(i)>stagnate_ratio*d_vec(max(i-1,1))
        u=u_old;
        B=B_func(u);
        d=B*u-c;
        idx_on=d>eps;
        fprintf('Penalty stagnated: penalty=%d in it=%d, norm(d)=%d, stagnated norm(d)=%d\n',rho,i,max(d(idx_on)),d_vec(i));
        break;
    end
    
    l_vypis = my_print(sprintf('It. %d: Penalty= %d, norm(d)=%d, geometry difference =%d \n',i,rho,max(d(idx_on)),res_geom),l_vypis );
    rho=rho*rho_step;
    u_old=u;
    B_old=B;
end
d_vec=d_vec(1:i);
if i==max_it
    fprintf('Maximum iterations reached: penalty=%d in it=%d, norm(d)=%d\n',rho,i,d_vec(end));
end
end

