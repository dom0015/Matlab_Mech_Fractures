FETI_problem_assembly

SMALSE_params.rel=1.0e-8;
SMALSE_params.rho0=1;
SMALSE_params.betarho=2;
SMALSE_params.Gama = 1;
SMALSE_params.M_start=1;
SMALSE_params.tol_to_update=1e-9;
SMALSE_params.maxiter_cg = 10000;
SMALSE_params.type='m';
SMALSE_params.print=true;

problem_setting.A=A;
problem_setting.A_plus=A_plus;
problem_setting.B_e=B_e;
problem_setting.b=b;
problem_setting.R=R;
problem_setting.c_e=c_e;
problem_setting.c_i=c_i;
problem_setting.B_iupdate=B_iupdate;
problem_setting.fracture_matrice=fracture_matrice;
problem_setting.plot_func2=plot_func2;

tic;
[D] = SMALSE_solver(problem_setting,SMALSE_params);
toc;