FETI_problem_assembly

SMALSE_params.rel=1.0e-8;
SMALSE_params.rho0=1;
SMALSE_params.betarho=2;
SMALSE_params.Gama = 1;
SMALSE_params.M_start=1;
SMALSE_params.tol_to_update=1e-9;
SMALSE_params.maxiter_cg = 10000;
SMALSE_params.type='m';
SMALSE_params.print=false;


frac_press={@(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)2*frac_press_val/5+0*x,@(x)2*frac_press_val/5+0*x
    @(x)frac_press_val/5.1+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x};

%% HYDRO
tocouple_aperture_independent
%% depends on fracture aperture
d = exp([-6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7])/100;
d = d(1:no_fractures);
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = d(i)*ones(lengths(i)-1,1);
end
figure(5); plot([mean(D{1}) mean(D{2}) mean(D{3}) mean(D{4}) mean(D{5}) mean(D{6}) mean(D{7}) mean(D{8}) ]); hold on
for i=1:20
    [PRESSURE,u0_]=tocouple_handle(D,no_fractures,mat_frac,fracture_matrice,...
        POINTS,intersections,alfa_inter,lengths,A,freeNode,b,u0);
    figure(6); plot([mean(PRESSURE{1}) mean(PRESSURE{2}) mean(PRESSURE{3}) mean(PRESSURE{4}) mean(PRESSURE{5}) mean(PRESSURE{6}) mean(PRESSURE{7}) mean(PRESSURE{8}) ]); hold on

    [problem_setting] = assembly_FETI_frac_rhs(problem_setting,PRESSURE);

    tic;
    [D] = SMALSE_solver(problem_setting,SMALSE_params);
    toc;
    
    for j=1:length(D)
        D{j}=D{j}+1e-5;
    end
    figure(5); plot([mean(D{1}) mean(D{2}) mean(D{3}) mean(D{4}) mean(D{5}) mean(D{6}) mean(D{7}) mean(D{8}) ]); hold on
end