function [Q,D,PRESSURE,ugrad,iter,x_elast,response_D] = coupled_solver_simple(params,hydro_problem,problem_setting,SMALSE_params)
global par_BiotWillis par_a0
hydro_problem.mat_frac=params;
no_fractures=hydro_problem.no_fractures;
lengths=hydro_problem.lengths;
d = 1e-4*ones(no_fractures,1);
% d = 0*ones(no_fractures,1);
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = d(i)*ones(lengths(i)-1,1);
end
% load('solution1_2e8.mat'); D{1}=solution1_2e8;
eps_coupling=SMALSE_params.eps_coupling;
res_d=[];
res_press=[];
alpha1=0.5;
alpha2=0.5;
flag=0;
tic;
% [PRESSURE__,~,ugrad,Q,PRESSURE]=coup.tocouple_handle_trhlina(D,hydro_problem);
[PRESSURE__,~,ugrad,Q,PRESSURE]=tocouple_handle_modif(D,hydro_problem);
D_old=D;
PRESSURE_old=PRESSURE;
ugrad_old=ugrad;
% ugrad_all=ugrad(:,1);
v_alpha=[];
kk=0;
rel_D_=[];
rel_P_=[];
response_D = cell2mat(D);
for i=1:no_fractures
    D{i}=0*D{i};
end
for i=1:SMALSE_params.coupling_iter
    %alpha=max(0.01,alpha*0.95);
%     [PRESSURE__,u0_,ugrad,Q,PRESSURE,frac_grad]=coup.tocouple_handle_trhlina(D,hydro_problem);
    D_=D;
    for j=1:no_fractures
        D_{j}=D{j}+par_a0;
    end
    [PRESSURE__,u0_,ugrad,Q,PRESSURE,frac_grad]=tocouple_handle_modif(D_,hydro_problem);
    tmp=cell2mat(PRESSURE__);
    res_press(:,i)=tmp(:);
%     ugrad_all(:,i)=ugrad(:,1);

    for j=1:size(PRESSURE,1)
        for kki=1:size(PRESSURE,2)
            PRESSURE{j,kki}=(PRESSURE{j,kki}*alpha1+(1-alpha1)*PRESSURE_old{j,kki});
        end
    end
    ugrad=ugrad*alpha1+(1-alpha1)*ugrad_old;
    
    [problem_setting] = feti.assembly_FETI_frac_rhs(problem_setting,PRESSURE,-ugrad*par_BiotWillis);
    [D,problem_setting,x_elast] = smalse.SMALSE_solver(problem_setting,SMALSE_params);
    for j=1:length(D)
        D{j}=max(D{j},0);%+1e-4;  
    end
    %temp = D{1}; temp(1)=temp(1)*2; temp(end)=temp(end)*2;
    response_D = [response_D cell2mat(D)];

    res_d(:,i)=cell2mat(D);
    for j=1:length(D)
        D{j}=(D{j}*alpha2+(1-alpha2)*D_old{j});
    end
%     for j=1:length(D)
%         temp=D{j}.^2; temp=log(temp);
%         temp0=D_old{j}.^2; temp0=log(temp0);
%         temp = temp*alpha2+(1-alpha2)*temp0;
%         temp = exp(temp); temp = sqrt(temp);
%         D{j}=temp;
%     end
    
    rel_P=sqrt(sum((res_press(:,1:end-1)-res_press(:,2:end)).^2,1))./sqrt(sum((res_press(:,2:end)).^2,1));
    rel_D=sqrt(sum((res_d(:,1:end-1)-res_d(:,2:end)).^2,1))./sqrt(sum((res_d(:,2:end)).^2,1));  
    
    % SMALSE precision to match coupled precision
    if ~isempty(rel_D)
        %SMALSE_params.rel=1e-9;%min(rel_D(end),rel_P(end))/100000;
    end
    % coupling parameters update
    
    if flag==1
        kk=kk+1;
        rel_D_=[rel_D_ rel_D(end)];
        rel_P_=[rel_P_ rel_P(end)];
    end
  
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
        fprintf('i=%d: alpha = %.4f , beta = %.2f || pres = %d || d = %d\n',i,alpha1,flag,rel_P(end),rel_D(end));
               
        if flag==1 && alpha1>=0.00001 && rel_P(end)<eps_coupling && rel_D(end)<eps_coupling
            break
        end
        if alpha1 <0.00001
            break
        end          
    end
    
    D_old=D;
    PRESSURE_old=PRESSURE;
    ugrad_old=ugrad;
    
    flag=1;
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
%     set(gca,'YScale','log')
    figure(6); title('alpha')
    hold on
    plot(v_alpha)
end
iter=length(rel_D_);
end

