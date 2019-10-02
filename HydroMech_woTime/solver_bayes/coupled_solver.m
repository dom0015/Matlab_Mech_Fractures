function [Q,D] = coupled_solver(params,hydro_problem,problem_setting,SMALSE_params)


SMALSE_params.rel=1.0e-2;
hydro_problem.mat_frac=params;
%% depends on fracture aperture
%d = exp([-6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7])/100;
no_fractures=hydro_problem.no_fractures;
lengths=hydro_problem.lengths;
d = 1e-4*ones(no_fractures,1);
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = d(i)*ones(lengths(i)-1,1);
end
D_old=D;
eps_coupling=SMALSE_params.eps_coupling;
res_d=[];
res_press=[];
alpha=0.2;
beta=1e-01;
tic;
for i=1:SMALSE_params.coupling_iter    
    [PRESSURE,~,ugrad,Q]=tocouple_handle(D,hydro_problem);
    res_press(:,i)=cell2mat(PRESSURE);
    
    if i>1
        for j=1:length(PRESSURE)
            PRESSURE{j}=(PRESSURE{j}*alpha+(1-alpha)*PRESSURE_old{j});
        end
        ugrad=ugrad*alpha+(1-alpha)*ugrad_old;
    end
    for j=1:length(PRESSURE)
        PRESSURE{j}=PRESSURE{j}*beta;
    end
    ugrad=beta*ugrad;
    [problem_setting] = assembly_FETI_frac_rhs(problem_setting,PRESSURE,-1*ugrad);
    
    [D,problem_setting] = SMALSE_solver(problem_setting,SMALSE_params);

    
    for j=1:length(D)
        D{j}=max(D{j},1e-10);
    end
    
    res_d(:,i)=cell2mat(D);
    if i>1
        for j=1:length(D)
            D{j}=(D{j}*alpha+(1-alpha)*D_old{j});
        end
    end
    
    tmp=sqrt(sum((res_press(:,1:end-1)-res_press(:,2:end)).^2,1))./sqrt(sum((res_press(:,2:end)).^2,1));
    tmp2=sqrt(sum((res_d(:,1:end-1)-res_d(:,2:end)).^2,1))./sqrt(sum((res_d(:,2:end)).^2,1));
    if i>1
        if SMALSE_params.print
        fprintf('i=%d: alpha = %.2f , beta = %.2f || pres = %d || d = %d\n',i,alpha,beta,tmp(end),tmp2(end));
        end
    end
    D_old=D;
    PRESSURE_old=PRESSURE;
    ugrad_old=ugrad;
    
    
    % SMALSE precision to match coupled precision
    if i>2
        SMALSE_params.rel=tmp2(end)/10;
    end
    % coupling parameters update
    alpha=min(0.75,alpha*1.1);
    beta=min(1,beta^(0.75)+1e-2);
    if i>1
        if beta==1 && alpha>=0.2 && tmp(end)<eps_coupling && tmp2(end)<eps_coupling
            break
        end
    end
end
end

