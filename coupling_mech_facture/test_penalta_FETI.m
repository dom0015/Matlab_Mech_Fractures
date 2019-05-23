
d_eps=1e-8;
d_min=0;
max_it=1000;
rho_step=2;
rho_step_limits=10;
stagnate_count=20;
alpha=0.75;


[problem_setting] = assembly_FETI_frac_rhs(problem_setting,PRESSURE,-1*ugrad);

A=problem_setting.A;
%problem_setting.A_plus=A_plus;
B_e=problem_setting.B_e;
b=problem_setting.b_frac;
%problem_setting.R=R;
c_e=problem_setting.c_e;
%problem_setting.c_i=c_i;
B_iupdate=problem_setting.B_iupdate;
fracture_matrice=problem_setting.fracture_matrice;
plot_func=problem_setting.plot_func2;
%problem_setting.mat_scale=mat_scale;






[u,i,d_vec,un_vec,geom_vec,rho,rhodtab] = penalta_FETI_noswitch(A,B_e,b,c_e,B_iupdate,d_min,max_it,d_eps,rho_step,rho_step_limits,stagnate_count,alpha);

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
B_i=B_iupdate(u);
[D_penalta] = construct_apertures(B_i*u,fracture_matrice);