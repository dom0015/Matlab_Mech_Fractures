
Nxy=101;
L1=1; L2=1;
sumbdomains_FETI=ceil(Nxy/10)^2;

mat_const=1e9;
frac_press_val=1;
frac_start_end={[0.2 0.4], [0.9 0.4];
                [0.2 0.2], [0.8 0.8]}; 
            
mat_omega_const=1e-15;
mat_frac_const=1e-6;
alfa_inter_const=1e-5;
                  
FETI_problem_assembly
tocouple_aperture_independent

SMALSE_params.rel=1.0e-2;
SMALSE_params.rho0=1;
SMALSE_params.betarho=2;
SMALSE_params.Gama = 1;
SMALSE_params.M_start=0.5;
SMALSE_params.tol_to_update=1e-4;
SMALSE_params.maxiter_cg = 100000;
SMALSE_params.type='m';
SMALSE_params.print=false;
problem_setting.B_i=[];
problem_setting.lambda_ker=[];

clearvars -except SMALSE_params problem_setting hydro_problem

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
eps_coupling=1e-5;
res_d=[];
res_press=[];
alpha=0.6;
beta=1e-1;
tic;
for i=1:200
    
    [PRESSURE,u0_,ugrad]=tocouple_handle(D,hydro_problem);
    res_press(:,i)=cell2mat(PRESSURE);
    if SMALSE_params.print
        figure(2);
        N=length(A);
        h = trisurf(ELEMENTS,POINTS(:,1),POINTS(:,2),u0_(1:N));
    end
    PRESSURE_=PRESSURE;
    ugrad_=ugrad;
    if i>1
        %alpha=max(0.1,min(1,alpha*norm(res_press(:,end))/norm(res_press(:,end)-res_press(:,end-1))/1e3));
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
    
    tic;
    [D,problem_setting] = SMALSE_solver(problem_setting,SMALSE_params);
    toc;
    
    for j=1:length(D)
        D{j}=max(D{j},1e-10);
    end
    D_=D;
    res_d(:,i)=cell2mat(D);
    if i>1
        for j=1:length(D)
            D{j}=(D{j}*alpha+(1-alpha)*D_old{j});
        end
    end
    
    figure(1)
    hold off
    plot(sqrt(sum((res_press(:,1:end-1)-res_press(:,2:end)).^2,1))./sqrt(sum((res_press(:,2:end)).^2,1)))
    tmp=sqrt(sum((res_press(:,1:end-1)-res_press(:,2:end)).^2,1))./sqrt(sum((res_press(:,2:end)).^2,1));
    hold on
    plot(sqrt(sum((res_d(:,1:end-1)-res_d(:,2:end)).^2,1))./sqrt(sum((res_d(:,2:end)).^2,1)))
    tmp2=sqrt(sum((res_d(:,1:end-1)-res_d(:,2:end)).^2,1))./sqrt(sum((res_d(:,2:end)).^2,1));
    set(gca,'YScale','log')
    drawnow
    D_old=D;
    PRESSURE_old=PRESSURE;
    ugrad_old=ugrad;
    fprintf('alpha=%d beta=%d\n',alpha,beta);
    if i>2
        SMALSE_params.rel=tmp2(end)/10;
    end
%     if i>2
%         if tmp(end)<=min(tmp)
%            alpha=min(1,alpha*1.1); 
%         else
%             alpha=max(0.1,alpha/1.1);
%         end
%     end
%     alpha
    %alpha=min(alpha*1.01,0.8);
%     if i>1
%         alpha=min(0.2,alpha^(0.99));
%     end
alpha=min(1,alpha*1.1);
    beta=min(1,beta^(0.5)+1e-1);
    if i>1
    if beta==1 && alpha>=0.2 && tmp(end)<eps_coupling && tmp2(end)<eps_coupling
        break
    end
    end
end
toc
[fig_id1,fig_id2,fig_id3] = plot_stresses(problem_setting.x_elast,problem_setting.NAPETI,problem_setting.sub_elem,problem_setting.sub_nodes,hydro_problem.POINTS,problem_setting.fracture_matrice);