
d_eps=1e-8;
d_min=0;
max_it=1000;
rho_step=2;
rho_step_limits=10;
stagnate_count=20;
alpha=0.75;

[u,i,d_vec,un_vec,geom_vec,rho] = penalta_FETI_noswitch(A,B_e,b,c_e,B_iupdate,d_min,max_it,d_eps,rho_step,rho_step_limits,stagnate_count,alpha);

figure
plot(d_vec,'r');
set(gca,'YScale','log')
hold on
plot(un_vec,'b');
plot(geom_vec,'g');

figure
plot(rho)
set(gca,'YScale','log')

plot_func(u,fracture_matrice);
