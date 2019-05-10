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


frac_press={@(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)2*frac_press_val/5+0*x,@(x)2*frac_press_val/5+0*x
    @(x)frac_press_val/5.1+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x};

[problem_setting] = assembly_FETI_frac_rhs(problem_setting,frac_press);

tic;
[D] = SMALSE_solver(problem_setting,SMALSE_params);
toc;