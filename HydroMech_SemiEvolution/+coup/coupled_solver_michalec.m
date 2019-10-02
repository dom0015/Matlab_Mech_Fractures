function [Q,D,PRESSURE,ugrad,iter,x_elast] = coupled_solver_michalec(params,hydro_problem,problem_setting,SMALSE_params)

SMALSE_params.rel=1.0e-16;
hydro_problem.mat_frac=params;
%% depends on fracture aperture
%d = exp([-6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7])/100;
no_fractures=hydro_problem.no_fractures;
lengths=hydro_problem.lengths;
d = 1e-3*ones(no_fractures,1);
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = d(i)*ones(lengths(i)-1,1);
end
eps_coupling=SMALSE_params.eps_coupling;
res_d=[];
res_press=[];
alpha=1;
beta=1e-2;
tic;
[PRESSURE__,~,ugrad,Q,PRESSURE]=coup.tocouple_handle_trhlina(D,hydro_problem);
D_old=D;
PRESSURE_old=PRESSURE;
ugrad_old=ugrad;
% ugrad_all=ugrad(:,1);
v_alpha=[];
kk=0;
rel_D_=[];
rel_P_=[];
alpha_all=[];
kkk=0;
for i=1:SMALSE_params.coupling_iter
    v_alpha=[v_alpha alpha];
    [PRESSURE__,u0_,ugrad,Q,PRESSURE,frac_grad]=coup.tocouple_handle_trhlina(D,hydro_problem);
    tmp=cell2mat(PRESSURE__);
    res_press(:,i)=tmp(:);
%     ugrad_all(:,i)=ugrad(:,1);

    for j=1:size(PRESSURE,1)
        for kki=1:size(PRESSURE,2)
            PRESSURE{j,kki}=(PRESSURE{j,kki}*alpha+(1-alpha)*PRESSURE_old{j,kki});
        end
    end
    ugrad=ugrad*alpha+(1-alpha)*ugrad_old;

%     for j=1:size(PRESSURE,1)
%         for kki=1:size(PRESSURE,2)
%             PRESSURE{j,kki}=PRESSURE{j,kki}*beta;
%         end
%     end
%     ugrad=beta*ugrad;
    
    [problem_setting] = feti.assembly_FETI_frac_rhs(problem_setting,PRESSURE,-ugrad);
    [D,problem_setting,x_elast] = smalse.SMALSE_solver(problem_setting,SMALSE_params);
%     [D,problem_setting,x_elast] = elasticity_solver(problem_setting,SMALSE_params);
    for j=1:length(D)
        D{j}=max(D{j},1e-5);  
    end
    
    
    
    res_d(:,i)=cell2mat(D);
    for j=1:length(D)
        D{j}=(D{j}*alpha+(1-alpha)*D_old{j});
    end
    
    
    rel_P=sqrt(sum((res_press(:,1:end-1)-res_press(:,2:end)).^2,1))./sqrt(sum((res_press(:,2:end)).^2,1));
    rel_D=sqrt(sum((res_d(:,1:end-1)-res_d(:,2:end)).^2,1))./sqrt(sum((res_d(:,2:end)).^2,1));  
    
    % SMALSE precision to match coupled precision
    if ~isempty(rel_D)
        SMALSE_params.rel=1e-16;%min(rel_D(end),rel_P(end))/100000;
    end
    % coupling parameters update
    
    if beta==1
        kk=kk+1;
        rel_D_=[rel_D_ rel_D(end)];
        rel_P_=[rel_P_ rel_P(end)];
        alpha_all=[alpha_all alpha];
        if kk>1
            ldD=(rel_D_(1)/rel_D_(end-1))^(1/(kk-1));
            ldP=(rel_P_(1)/rel_P_(end-1))^(1/(kk-1));
            cdD=(rel_D_(end-1)/rel_D_(end)).^(1/1);
            cdP=(rel_P_(end-1)/rel_P_(end)).^(1/1);
            
            if (ldD>cdD^(1.5))|| (ldP>cdP^(1.5))
                alpha=max(0,alpha*0.5);                
            end
            if (ldD<cdD)&&(ldP<cdP)
                alpha=min(1,alpha*1.1);
            end
%             tmp=min(1+cdD-ldD,1+cdP-ldP);
%             tmp=min(tmp,10);
%             tmp=max(tmp,0.1);
%             if tmp<1
%                 aa=0.8;
%                 alpha=alpha*1.1;%;aa*alpha+(1-aa)*alpha*tmp;
%                 alpha=min(alpha,1);
%             else
%                 aa=0.8;
%                 alpha=alpha*0.8;%aa*alpha+(1-aa)*alpha*tmp;
%                 alpha=max(alpha,0.001);
%             end
            
        end
    end
    %rel_D=rel_D./v_alpha(2:end);
    %rel_P=rel_P./v_alpha(2:end);
    beta=min(1,beta^(0.75)+1e-1);
  
    if SMALSE_params.print
        figure(13);
        N=length(hydro_problem.A);
        trisurf(hydro_problem.ELEMENTS,hydro_problem.POINTS(:,1),hydro_problem.POINTS(:,2),u0_(1:N),'LineStyle','none');
        figure(12); 
        plot(linspace(3,7,length(u0_(N+1:end))),u0_(N+1:end))
        grid on
        figure(11); %hold off
        %plot(linspace(3,7,length(-1e-5*(u0_(N+2:end)-u0_(N+1:end-1))/(0.1*sqrt(2)))),-1e-5*(u0_(N+2:end)-u0_(N+1:end-1))/(0.1*sqrt(2)))
        %hold on
        tmp=frac_grad{1};
        tmp(1)=tmp(1)/2;
        tmp(end)=tmp(end)/2;
        plot(linspace(3,7,length(-1e-5*tmp)),-1e-5*tmp)
    end
    
    if ~isempty(rel_D)
        fprintf('i=%d: alpha = %.4f , beta = %.2f || pres = %d || d = %d\n',i,alpha,beta,rel_P(end),rel_D(end));
               
        if beta==1 && alpha>=0.00001 && rel_P(end)<eps_coupling && rel_D(end)<eps_coupling
            break
        end
        if alpha <0.00001
            break
        end          
    end
    
    D_old=D;
    PRESSURE_old=PRESSURE;
    ugrad_old=ugrad;
    alpha=alpha*0.9999;
end

if SMALSE_params.print_couple
    figure(5)
    subplot(1,2,1)
    hold on
    plot(rel_P_,'k-','LineWidth',1)
    set(gca,'YScale','log')
    %figure(5)
    subplot(1,2,2)
    yyaxis right
    hold on
    plot(rel_D_,'k-','LineWidth',1)
    set(gca,'YScale','log')
    figure(6)
    hold on
    plot(v_alpha)
end
iter=length(rel_D_);
end

