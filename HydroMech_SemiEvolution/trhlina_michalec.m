global par_BiotWillis par_a0

%parametry Tom
par_permeabilita_trhliny = 1e-8;
par_dynamicka_viskozita = 0.001;
par_permeabilita_horniny = 1e-12*par_dynamicka_viskozita;
par_tloustka_trhliny = {1e-3,1e-3};
par_Lame_lambda = 1.8e9;
par_Lame_mu = 4.2e9;
par_BiotWillis = 0.89;
par_a0 = 1e-10;

% par_permeabilita_horniny = 1e-17;
% par_permeabilita_trhliny = 1e-8;
% par_dynamicka_viskozita = 0.001;
% par_tloustka_trhliny = 1e-3;
% par_Poissonovo_cislo = 0.3;
% par_Younguv_modul_horniny = 1e9;
% par_Lame_lambda = par_Younguv_modul_horniny*par_Poissonovo_cislo/(1+par_Poissonovo_cislo)/(1-2*par_Poissonovo_cislo);
% par_Lame_mu = par_Younguv_modul_horniny/2/(1+par_Poissonovo_cislo);
Nxy=101;
L1=10; L2=L1;
sumbdomains_FETI=100;

d_trhlina=par_tloustka_trhliny;
frac_press_val=1;
% frac_start_end={[0.3 0.3], [0.7 0.7]};
%frac_start_end={[0.1 0.5], [0.9 0.5]};
frac_start_end={[0.3 0.3], [0.7 0.7]
                [0.1 0.5], [0.9 0.5]};
fracture_direction={'d','h'};
no_fractures = 2;

mat_omega_const=par_permeabilita_horniny/par_dynamicka_viskozita;
mat_frac_const=par_permeabilita_trhliny/par_dynamicka_viskozita;
alfa_inter_const=1e-5;
cislo_ulohy=4;
n_windows=1;
DIRICHLET_PRESSURE=1e8;%5e8; % puvodne1e5; % 0,1,2,5

FETI_problem_assembly_michalec
hydro_preparation
hydro_problem.model=0;

SMALSE_params.rel=1.0e-9;
SMALSE_params.rho0=1;
SMALSE_params.betarho=2;
SMALSE_params.Gama = 1;
SMALSE_params.M_start=0.5;
SMALSE_params.tol_to_update=1e3;
SMALSE_params.maxiter_cg = 2000;
SMALSE_params.type='m';
SMALSE_params.print=false;
SMALSE_params.print_couple=true;
SMALSE_params.coupling_iter=100;
SMALSE_params.eps_coupling=1e-8;
problem_setting.B_i=[];
problem_setting.lambda_ker=[];
problem_setting.c_i=[];%d_trhlina;
% clearvars -except SMALSE_params problem_setting hydro_problem cislo_ulohy iters1 iters2 iters3 iters4 node_map_on
global S_iter
S_iter=[];
rng(1)
rand_add=randn(200,4);
iters=zeros(200,1);
for i=1:1
    params=hydro_problem.mat_frac;
    t=tic;
    [~,D,P,ug,iters(i),x_elast,response_D] = coup.coupled_solver_simple(params,hydro_problem,problem_setting,SMALSE_params);
    tt=tic;
    fprintf('Time of calculation = %.2f seconds.\n',double(tt-t)/1e6);
end

x_elast_demap=x_elast(node_map_on);
elast = [x_elast_demap(1:2:end), x_elast_demap(2:2:end)];
above = hydro_problem.fracture_matrice{1}.above_nodes;
above = [above(:,1); above(end,2)];
under = hydro_problem.fracture_matrice{1}.under_nodes;
under = [under(:,1); under(end,2)];
LEN = length(elast(above,1));
figure(103); plot(linspace(3,7,LEN),elast(above,1)); hold on; plot(linspace(3,7,LEN),elast(under,1))
grid on
%xlim([2.9, 7.1]); ylim([7.2, 8.7]*1e-5);
legend('up','down')
figure(104); plot(linspace(3,7,LEN),elast(above,2)); hold on; plot(linspace(3,7,LEN),elast(under,2))
grid on
%xlim([2.9, 7.1]); ylim([-18, -3]*1e-6);
legend('up','down')
NORM = -[1; -1]/sqrt(2);
figure(105); plot(linspace(3,7,LEN),elast(above,:)*NORM); hold on;
plot(linspace(3,7,LEN),elast(under,:)*NORM);
legend('up','down')
grid on
%xlim([2.9, 7.1]); ylim([-7.1, -5.6]*1e-5)
figure(13); colormap jet; axis equal; view(0,90); title('Hydraulic pressure')
figure(102); title('FETI subdomains')
figure(11); grid on; %xlim([2.9, 7.1]); ylim([0 0.44]*1e-6);
title('Darcy velocity in 1d fracture domain')
% figure; plot(linspace(3,7,LEN-1),D{1}); grid on; title("Aperture")
figure; imagesc(response_D); colorbar
figure; imagesc(log10(max(response_D,1e-5))); colorbar