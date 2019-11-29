function [Q,D,PRESSURE,ugrad,iter,x_elast,response_D] = coupled_solver_adaptive1(hydro_problem,elast_problem,SMALSE_params,initial_aperture)
par_BiotWillis = hydro_problem.par_BiotWillis;
no_fractures=hydro_problem.no_fractures;
lengths=hydro_problem.lengths;
eps_coupling=SMALSE_params.eps_coupling;

D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = initial_aperture(i)/10*ones(lengths(i)-1,1);
end

[~,u_old_single,ugrad,Q,PRESSURE]=coup.hydro_basicStationary(D,hydro_problem);
u_old={u_old_single};
ugrad_all={ugrad};
PRESSURE_all={PRESSURE};
[elast_problem] = feti.assembly_FETI_frac_rhs(elast_problem,PRESSURE,-ugrad*par_BiotWillis);
[D,elast_problem,x_elast] = smalse.SMALSE_solver(elast_problem,SMALSE_params);
D_all={D};

shifts=[0];
cas_krok=[0];
hydro_change=[0];
alphas=[0];
trasholds=[0];
smooth_distances=[0];
D_log_distance_all=[0];

response_D = cell2mat(D);
response_D_smooth=cell2mat(D);
trashold=1;
trashold_D=0.5;
alpha=1;
beta1=10;
beta2=2;
load('D_fin2.mat')
i=1;
iter=0;
tmp_prec=1;
while iter<=SMALSE_params.coupling_iter
    i=i+1;
    iter=iter+1;
    hydro_problem.const_delta_t=1e-1;
    [~,u0_,ugrad,Q,PRESSURE,~]=coup.hydro_basicSemiEvolution(D_all{i-1},hydro_problem,u_old{i-1});
    shift=hydro_shift(ugrad_all{i-1},ugrad);
    shifts(i)=shift;
    
    cas_krok(i)=hydro_problem.const_delta_t;
    fprintf('  IT: %d ,Shift: %d, delta T: %d, trasholdD: %d \n',i,shift,hydro_problem.const_delta_t,trashold_D)
    
    u_old{i} = u0_;
    ugrad_all{i}=ugrad;
    PRESSURE_all{i}=PRESSURE;
    
    [elast_problem] = feti.assembly_FETI_frac_rhs(elast_problem,PRESSURE_all{i},-ugrad_all{i}*par_BiotWillis);
    [D,elast_problem,x_elast,ncg] = smalse.SMALSE_solver(elast_problem,SMALSE_params);
    response_D(:,i) =cell2mat(D);
    
    if ncg<50
        SMALSE_params.rel=max(SMALSE_params.rel*0.5,1e-12);
    end
    
    [~,~,ugrad_stationary,~,~]=coup.hydro_basicStationary(D,hydro_problem);
    sta_dist=hydro_shift(ugrad_all{i},ugrad_stationary);
    hydro_change(i)=sta_dist;
    fprintf('  Stationary distance: %d',sta_dist)
    
    if iter>2
        alpha=min(trashold_D/mech_shift(response_D_smooth(:,i-1),response_D(:,i)),1);
    end

    
    
    D=smooth_D(D_all{i-1},D,alpha);
    D_all{i}=D;
    response_D_smooth(:,i)=cell2mat(D);
    
    smooth_distance=mech_shift(response_D_smooth(:,i-1),response_D_smooth(:,i));
    smooth_distances(i)=smooth_distance;
    
    
    D_log_distance_loc=mech_shift(response_D(:,i-1),response_D(:,i));
    D_log_distance_all(i)=D_log_distance_loc;
    
    fprintf('--- D_now: %d, alpha: %d beta1: %d beta2 %d',D_log_distance_loc,alpha,beta1,beta2)
    
    alphas(i)=alpha;
    
    fprintf('\n')
    if sta_dist/min(hydro_problem.const_delta_t,1)<eps_coupling && (D_log_distance_loc)<0.1
        break
    end
    
    trashold_Ds(i)=trashold_D;
    
    
%     trashold_D=trashold_D/1.05;
%     hydro_problem.const_delta_t=hydro_problem.const_delta_t*2;
%     if alpha>0.2
%         hydro_problem.const_delta_t=hydro_problem.const_delta_t*2;
%     end
%     if alpha<0.1
%         if sta_dist< D_log_distance_loc
%             hydro_problem.const_delta_t=hydro_problem.const_delta_t/2;
%         end
%         trashold_D=trashold_D*1.2;
%     end
end

fprintf('iterations: %d\n',i);
% figure
% plot(time_change)
% figure
% plot(D_all_dist)
% hold on
% plot(D_all_old_dist)

figure
plot(hydro_change)
hold on
plot(D_log_distance_all)
hold on
plot(1./cas_krok)
plot(alphas)
plot(smooth_distances)
plot(trashold_Ds)
set(gca,'YScale','log')


figure
plot(shifts)
hold on
plot(smooth_distances)
set(gca,'YScale','log')

figure
imagesc(response_D)
iter=i;

figure
trisurf(hydro_problem.ELEMENTS,hydro_problem.POINTS(:,1),hydro_problem.POINTS(:,2),u0_(1:length(hydro_problem.POINTS)),'LineStyle','none');
end

function shift=hydro_shift(u_old,u_new)
shift=sort(abs(u_old(:)-u_new(:)))/1e6;%/max(abs(u_new(:)));
shift=shift(ceil(end*0.95));
end

function shift=mech_shift(u_old,u_new)
shift=sort(abs(log(u_new)-log(u_old)));
shift=shift(ceil(end*0.95));
end

function D=smooth_D(D_old,D,alpha)
for j=1:length(D)
    D{j}=exp((log(max(D{j},1e-12))*alpha+(1-alpha)*log(max(D_old{j},1e-12))));
end
end