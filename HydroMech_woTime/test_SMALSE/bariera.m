function [u,i] = bariera(A,b,B_func,d_min,max_it,rho_step,d_eps,pl,u_0)

B=-B_func(u_0);
B_t=B;
idx_on=B*u_0<1e-4;
B=B(idx_on,:);
M=[A B';B sparse(size(B,1),size(B,1))];
f=[b;1e-3*ones(size(B,1),1)];
B=B_t;
tmp0=M\f;
u=tmp0(1:length(b));
pl(u);
rho0=1e-4;


H_P=B'*diag(rho0*sparse((B*u).^(-2)))*B;
grad_P=-B'*(rho0./(B*u));

c=d_min*ones(size(B,1),1);
u_old=u;
l_vypis=0;
for i=1:max_it
    alpha=0.5;
    beta=1;
    u_n=((A+H_P)\(A*u+grad_P-b));
    u_t=u-alpha*u_n;
    while sum(B*u_t-c<rho0)>0 || max(abs(alpha*u_n))>0.05
        alpha=alpha/2;
        u_t=u-alpha*u_n;
    end
    u=u_t;
    
    B_new=-B_func(u);
    B_t=(1-beta)*B+beta*B_new;
    
    while sum(B_t*u-c<rho0)>0
        if beta<0.2
            [k] = proj_B_larger(u,B_t,c,rho0);
            u=u+k;
            pl(u);
            break
        end
        beta=beta/2;
        B_t=(1-beta)*B+beta*B_new;

    end
    B=B_t;
    d=B*u-c;
    change=norm(u-u_old);
    if change<(rho0*10) && beta > 0.1
        rho0=rho0/rho_step;
        if min(d)<d_eps 
            break
        end
    end
    rho=rho0;
    tmp1=max(abs(rho*((B*u).^(-2))));
    tmp2=max(abs(rho./(B*u)));
    H_P=B'*sparse(rho*diag((B*u).^(-2)))*B;
    grad_P=-B'*(rho./(B*u));
    
    u_old=u;
    l_vypis = my_print(sprintf('It. %d: Bariera= %d, norm(d)=%d, norm(u)=%d, beta=%d, H=%d, g=%d\n',i,rho,min(d),change,beta,tmp1,tmp2),l_vypis );
end
end