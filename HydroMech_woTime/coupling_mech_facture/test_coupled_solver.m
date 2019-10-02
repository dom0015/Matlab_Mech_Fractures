
Nxy=101;
L1=1; L2=1;
sumbdomains_FETI=ceil(Nxy/10)^2;

mat_const=1e9;
frac_press_val=1;
frac_start_end={[0.3 0.2], [0.9 0.2]%
     [0.3 0.8], [0.9 0.8]%
     [0.1 0.3], [0.7 0.9]%
     [0.4 0.1], [0.4 0.7]%
     };
            
mat_omega_const=1e-15;
mat_frac_const=1e-6;
alfa_inter_const=1e-5;
cislo_ulohy=4;
n_windows=32;
                  
FETI_problem_assembly
tocouple_aperture_independent

SMALSE_params.rel=1.0e-2;
SMALSE_params.rho0=1;
SMALSE_params.betarho=2;
SMALSE_params.Gama = 1;
SMALSE_params.M_start=0.5;
SMALSE_params.tol_to_update=1e-4;
SMALSE_params.maxiter_cg = 2000;
SMALSE_params.type='m';
SMALSE_params.print=false;
SMALSE_params.coupling_iter=41; 
SMALSE_params.eps_coupling=1e-6;
problem_setting.B_i=[];
problem_setting.lambda_ker=[];

clearvars -except SMALSE_params problem_setting hydro_problem cislo_ulohy

params=[1 1 1 1]*10^(-6);

t=tic;
[Q,D] = coupled_solver(params,hydro_problem,problem_setting,SMALSE_params);
tt=tic;
fprintf('Time of calculation = %.2f seconds.\n',double(tt-t)/1e6);

figure(11)
subplot(1,2,1)
hold on
plot(cell2mat(D));
subplot(1,2,2)
hold on
plot(Q)
max(Q)
save(['uloha' num2str(cislo_ulohy) '.mat'],'problem_setting','hydro_problem','SMALSE_params','D','Q')
