
%% MESH PARAMETERS
shared_data.Nxy=101;
shared_data.L1=10; shared_data.L2=10;
shared_data.frac_start_end={[0.3 0.3], [0.7 0.7]
                [0.1 0.5], [0.9 0.5]};
            
%% ELAST PARAMETERS
elast_problem.sumbdomains_FETI=100;
elast_problem.par_tloustka_trhliny = 1e-3;
elast_problem.par_Lame_lambda = 1.8e9;
elast_problem.par_Lame_mu = 4.2e9;

%% HYDRO PARAMETERS
hydro_problem.cislo_ulohy=4;
hydro_problem.n_windows=1;
hydro_problem.DIRICHLET_PRESSURE=1e8; % puvodne1e5; % 
hydro_problem.hydro_model=0;
hydro_problem.alfa_inter_const=1e-5;
par_permeabilita_trhliny = 1e-8;
par_dynamicka_viskozita = 0.001;
par_permeabilita_horniny = 1e-15;%permeabilita
hydro_problem.mat_omega_const=par_permeabilita_horniny/par_dynamicka_viskozita;%hydraulicka konduktivita
hydro_problem.mat_frac_const=par_permeabilita_trhliny/par_dynamicka_viskozita;

%% OTHER INPUT PARAMETERS
par_storativity = 0.1449e-9;
hydro_problem.const_cs_domain = hydro_problem.mat_omega_const;
hydro_problem.const_cs_fracture = hydro_problem.mat_frac_const*elast_problem.par_tloustka_trhliny^2;
hydro_problem.par_BiotWillis = 0.89;
hydro_problem.par_a0 = 1e-6;
hydro_problem.const_delta_t = 1;
%% SOLVER PARAMETERS
SMALSE_params.rel=1.0e-4;
SMALSE_params.rho0=1;
SMALSE_params.betarho=2;
SMALSE_params.Gama = 1;
SMALSE_params.M_start=0.5;
SMALSE_params.tol_to_update=1e3;
SMALSE_params.maxiter_cg = 3000;
SMALSE_params.type='m';
SMALSE_params.print=false;
SMALSE_params.print_couple=true;
SMALSE_params.coupling_iter=70;
SMALSE_params.eps_coupling=1e-3;

[elast_problem,shared_data] = elast_preparation(elast_problem,shared_data);
hydro_problem = hydro_preparation( hydro_problem,shared_data );
initial_aperture = 1e-4*ones(hydro_problem.no_fractures,1);
[Q,D,PRESSURE,ugrad,iter,x_elast,response_D] = coup.coupled_solver_adaptiveSemiEvolution(hydro_problem,elast_problem,SMALSE_params,initial_aperture);

%viz.fracture_displacement(elast_problem,x_elast)

figure; imagesc(response_D); colorbar; title("All apertures (cols) in coupling iterations")
figure; imagesc(log10(max(response_D,1e-5))); colorbar; title("Log of all apertures (cols) in coupling iterations")