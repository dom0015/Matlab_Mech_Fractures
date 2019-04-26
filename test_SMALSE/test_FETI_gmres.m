
r=1e-3;


rho=1e-6;
n1=size(B_e,1);
n2=size(B_i,1);
O_e=sparse(n1,n1);
I_i=1/rho^2*speye(n2,n2);
O_ie=sparse(n1,n2);

D=[O_e     O_ie;
    O_ie'    I_i];
prec=@(x)preconditioner_FETI(x,A,A_plus,B,r);

A_full=[A       B_e'   B_i';...
        B_e    O_e     O_ie;
        B_i O_ie'    -I_i];


b_full=[b;c_e;c_i];

[x,flag,relres,iter,resvec]=gmres(A_full,b_full,100,1e-12,1,prec);
plot(resvec)
set(gca,'YScale','log')

plot_func(x,fracture_matrice);