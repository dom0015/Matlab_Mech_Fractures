function [Q,D,PRESSURE,ugrad,iter,x_elast,response_D] = coupled_solver_adaptive(hydro_problem,elast_problem,SMALSE_params,initial_aperture)
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
trashold_D=1;
alpha=1;

i=1;
iter=0;
tmp_prec=1;
while iter<=SMALSE_params.coupling_iter
    i=i+1;
    iter=iter+1;
    
    [~,u0_,ugrad,Q,PRESSURE,~]=coup.hydro_basicSemiEvolution(D_all{i-1},hydro_problem,u_old{i-1});
    shift=hydro_shift(ugrad_all{i-1},ugrad);
    shifts(i)=shift;
    
    if shift<trashold*0.9
        hydro_problem.const_delta_t=min(hydro_problem.const_delta_t*2,1/max(trashold,trashold_D));
    end
    
    while shift>=trashold
        hydro_problem.const_delta_t=hydro_problem.const_delta_t/2;
        [~,u0_,ugrad,Q,PRESSURE,~]=coup.hydro_basicSemiEvolution(D_all{i-1},hydro_problem,u_old{i-1});
        shift=hydro_shift(ugrad_all{i-1},ugrad);
        fprintf('    Correction: Shift: %d, delta T: %d \n',shift,hydro_problem.const_delta_t)
    end
    
    cas_krok(i)=hydro_problem.const_delta_t;
    fprintf('  IT: %d ,Shift: %d, delta T: %d, thashold: %d, trasholdD: %d \n',i,shift,hydro_problem.const_delta_t,trashold,trashold_D)
    
    u_old{i} = u0_;
    ugrad_all{i}=ugrad;
    PRESSURE_all{i}=PRESSURE;
    
    
    [elast_problem] = feti.assembly_FETI_frac_rhs(elast_problem,PRESSURE_all{i},-ugrad_all{i}*par_BiotWillis);
    [D,elast_problem,x_elast,ncg] = smalse.SMALSE_solver(elast_problem,SMALSE_params);
    response_D(:,i) =cell2mat(D);
    
    if ncg<100
        SMALSE_params.rel=max(SMALSE_params.rel*0.5,1e-12);
    end
    
    [~,~,ugrad_stationary,~,~]=coup.hydro_basicStationary(D,hydro_problem);
    sta_dist=hydro_shift(ugrad_all{i},ugrad_stationary);
    hydro_change(i)=sta_dist;
    fprintf('  Stationary distance: %d',sta_dist)
    
    if iter>2
        %alpha_tmp=max(smooth_distances(i-1),sta_dist/2)/mech_shift(response_D_smooth(:,i-1),response_D(:,i));
        %alpha_tmp=trashold/mech_shift(response_D_smooth(:,i-1),response_D(:,i));
        %alpha=max(min(alpha_tmp,0.1),0.1);
        alpha=min(trashold_D/mech_shift(response_D_smooth(:,i-1),response_D(:,i)),0.5);
        %alpha=0.3;
    end
    
%     if alpha>0.1
%         trashold_D=trashold_D*0.75;
%     end
%     if alpha<0.075
%         hydro_problem.const_delta_t=hydro_problem.const_delta_t/10;
%     end
    
    
    D=smooth_D(D_all{i-1},D,alpha);
    D_all{i}=D;
    response_D_smooth(:,i)=cell2mat(D);
    
    smooth_distance=mech_shift(response_D_smooth(:,i-1),response_D_smooth(:,i));
    smooth_distances(i)=smooth_distance;
    
    
    D_log_distance_loc=mech_shift(response_D(:,i-1),response_D(:,i));
    D_log_distance_all(i)=D_log_distance_loc;
    
    fprintf('--- D_now: %d, alpha: %d',D_log_distance_loc,alpha)
    
    
    alphas(i)=alpha;
    trasholds(i)=trashold;
    %trashold=max(exp(log(trashold)*0.5+(log(sta_dist))*0.5),1e-6);
    
    
    if iter>1
        if D_log_distance_loc>mech_shift(response_D(:,i-2),response_D(:,i))
            trashold_D=trashold_D*0.5;
        else
            trashold_D=max(exp(log(trashold_D)*0.75+log(min(D_log_distance_loc/alpha,1)/2)*0.25),1e-6);
        end
    end
   % trashold_D=max(exp(log(trashold_D)*0.75+log(min(smooth_distance/alpha,1)/2)*0.25),1e-6);
    
%     if iter>3
%         if D_log_distance_all(i)>trashold
% %            trashold=trashold*0.5;
% %            alpha=alpha*0.8;
%             hydro_problem.const_delta_t=hydro_problem.const_delta_t/2;
% %              D_all{i}=D_all{i-1};
% %              D_log_distance_all(i)=D_log_distance_all(i-1);
% %              ugrad_all{i}=ugrad_all{i-1};
% %              u_old{i}=u_old{i-1};
% %              response_D(:,i)=response_D(:,i-1);
%             fprintf('XXXX')
%         else
%            hydro_problem.const_delta_t=min(hydro_problem.const_delta_t*10,1e5);
%         end
%         if D_log_distance_all(i)<min(D_log_distance_all(1:i-1))
%             %trashold=trashold*1.1;
%             %hydro_problem.const_delta_t=hydro_problem.const_delta_t*2;
%             %alpha=min(alpha*1.05,0.3);
%         end
%     end
    
    fprintf('\n')
    if sta_dist/min(hydro_problem.const_delta_t,1)<eps_coupling && (D_log_distance_loc/alpha)<eps_coupling*3
        break
    end
    
%     if iter>10
%     if D_log_distance_all(i)>min(0.9*D_log_distance_all(i-1),0.5)
%         hydro_problem.const_delta_t=hydro_problem.const_delta_t/10;
%     else
%         hydro_problem.const_delta_t=hydro_problem.const_delta_t*1.5;
%     end
%     
%     else
%         
%     hydro_problem.const_delta_t=hydro_problem.const_delta_t*1.5;
%     end
%     hydro_problem.const_delta_t=hydro_problem.const_delta_t*2;
    
%     if  sta_dist<trashold 
%         trashold=trashold/1.1
%     end
trashold=max(exp(log(trashold)*0.75+(log(sta_dist/2))*0.25),1e-6);
trashold_Ds(i)=trashold_D;
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
plot(trasholds)
plot(smooth_distances)
plot(trashold_Ds)
set(gca,'YScale','log')

figure
imagesc(response_D)
iter=i;
end

function shift=hydro_shift(u_old,u_new)
shift=sort(abs(u_old(:)-u_new(:)))/max(abs(u_new(:)));
shift=shift(ceil(end*0.95));
end

function shift=mech_shift(u_old,u_new)
shift=sort(abs(log(u_new)-log(u_old)));
shift=shift(ceil(end*0.95));
end

function D=smooth_D(D_old,D,alpha)
for j=1:length(D)
    D{j}=exp((log(max(D{j},1e-6))*alpha+(1-alpha)*log(max(D_old{j},1e-6))));
end
end