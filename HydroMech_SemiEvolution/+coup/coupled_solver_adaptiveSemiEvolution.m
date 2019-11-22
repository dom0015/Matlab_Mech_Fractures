function [Q,D,PRESSURE,ugrad,iter,x_elast,response_D] = coupled_solver_SemiEvolution(hydro_problem,elast_problem,SMALSE_params,initial_aperture)
par_BiotWillis = hydro_problem.par_BiotWillis;
par_a0 = hydro_problem.par_a0;
no_fractures=hydro_problem.no_fractures;
lengths=hydro_problem.lengths;
eps_coupling=SMALSE_params.eps_coupling;

if nargin<4
    initial_aperture = 1e-4*ones(no_fractures,1);
end
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = initial_aperture(i)*ones(lengths(i)-1,1);
end
%load('D_fin1.mat')
res_press=[];
tic;
[PRESSURE__,u_old,ugrad,Q,PRESSURE]=coup.hydro_basicStationary(D,hydro_problem);
[elast_problem] = feti.assembly_FETI_frac_rhs(elast_problem,PRESSURE,-ugrad*par_BiotWillis);
[D,elast_problem,x_elast] = smalse.SMALSE_solver(elast_problem,SMALSE_params);
for j=1:length(D)
    D{j}=max(D{j},par_a0);
end

D_old=D;
% elast_div=ugrad(:,1)*0; elast_div_diff=elast_div;
kk=0;
rel_D_=[];
rel_P_=[];
kkk=[];
response_D = cell2mat(D);
response_D_smooth=[];
%u_old = 0*u_old;
trashold=1;
alpha=0.3;
max_alpha=0.3;
for i=1:SMALSE_params.coupling_iter
    ugrad_old=ugrad;
    [PRESSURE__,u0_,ugrad,Q,PRESSURE,frac_grad]=coup.hydro_basicSemiEvolution...
        (D,hydro_problem,u_old);%,elast_div_diff);
    
    shift=max(abs(ugrad_old(:)-ugrad(:)))/max(abs(ugrad(:)));
    kkk=[kkk shift];
    
    
    if shift<trashold/2 
        hydro_problem.const_delta_t=min(hydro_problem.const_delta_t*2,1e6);
    end
    cas_krok(i)=hydro_problem.const_delta_t;
    fprintf('  Shift: %d, delta T: %d, thashold: %d\n',shift,hydro_problem.const_delta_t,trashold)
    tmp_time=hydro_problem.const_delta_t;
    while shift>=trashold %& hydro_problem.const_delta_t>1e-3 %& hydro_problem.const_delta_t > tmp_time/10
        hydro_problem.const_delta_t=hydro_problem.const_delta_t/3;
        [PRESSURE__,u0_,ugrad,Q,PRESSURE,frac_grad]=coup.hydro_basicSemiEvolution(D,hydro_problem,u_old);
        shift=max(abs(ugrad_old(:)-ugrad(:)))/max(abs(ugrad(:)));
        fprintf('    Correction: Shift: %d, delta T: %d \n',shift,hydro_problem.const_delta_t)
    end
    
    u_old = u0_;
    tmp=cell2mat(PRESSURE__);
    res_press(:,i)=tmp(:);

    [elast_problem] = feti.assembly_FETI_frac_rhs(elast_problem,PRESSURE,-ugrad*par_BiotWillis);
    [D,elast_problem,x_elast,ncg] = smalse.SMALSE_solver(elast_problem,SMALSE_params);
    if ncg<500
        SMALSE_params.rel=SMALSE_params.rel*0.5;
    end
    
    [~,~,ugrad_stationary,~,~]=coup.hydro_basicStationary(D,hydro_problem);    
    sta_dist=max(abs(ugrad_stationary(:)-ugrad(:)))/max(abs(ugrad(:)));
    fprintf('  Stationary distance: %d',sta_dist)
    
    time_change(i)=sta_dist/hydro_problem.const_delta_t;
    hydro_change(i)=sta_dist;
    response_D = [response_D cell2mat(D)];
    if i>1
        tmp1=log(max(cell2mat(D),par_a0));
        tmp2=log(cell2mat(D_old));
        tmp3=mean(abs(tmp1-tmp2));
        tmp3all(i)=tmp3;
        alpha_tmp=smooth_distance/tmp3;
        %alpha=min(alpha,alpha_tmp);
        alpha=max(min(alpha_tmp*2,max_alpha),0.000000001);
        %alpha=1;%max_alpha;
    end
    
    for j=1:length(D)
        D{j}=exp((log(max(D{j},par_a0))*alpha+(1-alpha)*log(D_old{j})));
    end
    response_D_smooth = [response_D_smooth cell2mat(D)];
    smooth_distance=mean(abs(log(cell2mat(D))-log(cell2mat(D_old))));
    smooth_distances(i)=smooth_distance;
    D_old=D;
    D_log_distance_loc=max(abs(log(response_D(:,end))-log(response_D(:,end-1))));
    if i>1
        D_log_distance_old=min(max(abs(log(response_D(:,end))-log(response_D(:,max(1,end-5):end-2)))));
        D_all_dist(i)=D_log_distance_loc;
        D_all_old_dist(i)=D_log_distance_old;
        if D_log_distance_loc>D_log_distance_old
            max_alpha=max_alpha*0.5;
            trashold=trashold*0.75;
        end
        if D_log_distance_loc<D_log_distance_old
            max_alpha=max_alpha/0.75;
        end
        
        max_alpha=max(min(max_alpha,0.3),0.00001);
        fprintf('--- D_now: %d D_old: %d, alpha: %d max=%d',D_log_distance_loc,D_log_distance_old,alpha,max_alpha)
    end
    fprintf('\n')
    if hydro_problem.const_delta_t>1 && sta_dist<SMALSE_params.eps_coupling && D_log_distance_loc/alpha<SMALSE_params.eps_coupling
        break
    end
    alphas(i)=alpha;
    trasholds(i)=trashold;
    trashold=max(exp(log(trashold)*0.75+(log(sta_dist))*0.25),1e-6);
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
plot(D_all_dist)
hold on
plot(1./cas_krok)
plot(alphas)
plot(trasholds)
plot(smooth_distances)

figure 
imagesc(response_D_smooth)
iter=i;
end