

SMALSE_params.rel=1.0e-12;
SMALSE_params.rho0=1;
SMALSE_params.betarho=2;
SMALSE_params.Gama = 1;
SMALSE_params.M_start=0.5;
SMALSE_params.tol_to_update=1e-8;
SMALSE_params.maxiter_cg = 100000;
SMALSE_params.type='m';
SMALSE_params.print=true;
problem_setting.B_i=[];
problem_setting.lambda_ker=[];




[problem_setting] = assembly_FETI_frac_rhs(problem_setting,PRESSURE,-1*ugrad);
[D,problem_setting] = SMALSE_solver(problem_setting,SMALSE_params);