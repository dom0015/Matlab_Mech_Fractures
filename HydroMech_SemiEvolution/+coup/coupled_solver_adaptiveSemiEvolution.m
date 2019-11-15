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
%load('D_fin.mat')
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
        hydro_problem.const_delta_t=min(hydro_problem.const_delta_t*2,1e10);
    end
    cas_krok(i)=hydro_problem.const_delta_t;
    fprintf('  Shift: %d, delta T: %d, thashold: %d\n',shift,hydro_problem.const_delta_t,trashold)

    while shift>=trashold
        hydro_problem.const_delta_t=hydro_problem.const_delta_t*0.5;
        [PRESSURE__,u0_,ugrad,Q,PRESSURE,frac_grad]=coup.hydro_basicSemiEvolution(D,hydro_problem,u_old);
        shift=max(abs(ugrad_old(:)-ugrad(:)))/max(abs(ugrad(:)));
        fprintf('    Correction: Shift: %d, delta T: %d \n',shift,hydro_problem.const_delta_t)
    end
      
    u_old = u0_;
    tmp=cell2mat(PRESSURE__);
    res_press(:,i)=tmp(:);
    
    [elast_problem] = feti.assembly_FETI_frac_rhs(elast_problem,PRESSURE,-ugrad*par_BiotWillis);
    [D,elast_problem,x_elast,ncg] = smalse.SMALSE_solver(elast_problem,SMALSE_params);
    if ncg<100
        SMALSE_params.rel=SMALSE_params.rel*0.5;
    end
    
    [~,~,ugrad_stationary,~,~]=coup.hydro_basicStationary(D,hydro_problem);
    sta_dist=max(abs(ugrad_stationary(:)-ugrad(:)))/max(abs(ugrad(:)));
    fprintf('  Stationary distance: %d',sta_dist)
    trashold=max(exp(log(trashold)*0.5+(log(sta_dist)-log(1))*0.5),1e-3);
    time_change(i)=sta_dist/hydro_problem.const_delta_t;
    hydro_change(i)=sta_dist;
    response_D = [response_D cell2mat(D)];
    if i>1
        tmp1=log(cell2mat(D));
        tmp2=log(cell2mat(D_old));
        tmp3=max(abs(tmp1-tmp2));
        alpha_tmp=smooth_distance/tmp3;
        %alpha=min(alpha,alpha_tmp);
        alpha=max(min(alpha_tmp,max_alpha),0.05);
    end
    
    for j=1:length(D)
        D{j}=exp((log(max(D{j},par_a0))*alpha+(1-alpha)*log(D_old{j})));
    end
    smooth_distance=max(abs(log(cell2mat(D))-log(cell2mat(D_old))));
    D_old=D;
    D_log_distance_loc=max(abs(log(response_D(:,end))-log(response_D(:,end-1))));
    if i>1
        D_log_distance_old=max(abs(log(response_D(:,end))-log(response_D(:,end-2))));
        D_all_dist(i)=D_log_distance_loc;
        D_all_old_dist(i)=D_log_distance_old;
        if D_log_distance_loc>D_log_distance_old
            max_alpha=max_alpha*0.8;
        end
        if D_log_distance_loc*1.2<D_log_distance_old
            max_alpha=max_alpha*1.05;
        end
        
       max_alpha=max(min(max_alpha,0.3),0.01);
        
        if D_log_distance_loc>1.2*D_log_distance_old
            trashold=trashold/2;
        end
        fprintf('--- D_now: %d D_old: %d, alpha: %d max=%d',D_log_distance_loc,D_log_distance_old,alpha,max_alpha)
    end
    fprintf('\n')
    if hydro_problem.const_delta_t>1 && sta_dist<SMALSE_params.eps_coupling && D_log_distance_loc/alpha<SMALSE_params.eps_coupling
        break
    end
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



iter=i;
end