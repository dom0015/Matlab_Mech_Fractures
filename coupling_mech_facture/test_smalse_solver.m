FETI_problem_assembly

SMALSE_params.rel=1.0e-9;
SMALSE_params.rho0=1;
SMALSE_params.betarho=2;
SMALSE_params.Gama = 1;
SMALSE_params.M_start=0.5;
SMALSE_params.tol_to_update=1e-12;
SMALSE_params.maxiter_cg = 100000;
SMALSE_params.type='m';
SMALSE_params.print=true;
problem_setting.B_i=[];
problem_setting.lambda_ker=[];

frac_press={@(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)2*frac_press_val/5+0*x,@(x)2*frac_press_val/5+0*x
    @(x)frac_press_val/5.1+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x
    @(x)frac_press_val/5+0*x,@(x)frac_press_val/5+0*x};

%% HYDRO
tocouple_aperture_independent
%% depends on fracture aperture
%d = exp([-6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7])/100;
d = 1e-1*ones(no_fractures,1);
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = d(i)*ones(lengths(i)-1,1);
end
D_old=D;

res_d=[];
res_press=[];
alpha=0.4;
for i=1:50

    [PRESSURE,u0_]=tocouple_handle(D,no_fractures,mat_frac,fracture_matrice,...
        POINTS,intersections,alfa_inter,lengths,A,freeNode,b,u0);
    res_press(:,i)=cell2mat(PRESSURE);
    if SMALSE_params.print
        figure(2);
        N=length(A);
        h = trisurf(ELEMENTS,POINTS(:,1),POINTS(:,2),u0_(1:N));
    end
    if i>1
        for j=1:length(PRESSURE)
            PRESSURE{j}=(PRESSURE{j}*alpha+(1-alpha)*PRESSURE_old{j});
        end  
    end
    [problem_setting] = assembly_FETI_frac_rhs(problem_setting,PRESSURE);

    tic;
    [D,problem_setting] = SMALSE_solver(problem_setting,SMALSE_params);
    toc;
    
    for j=1:length(D)
        D{j}=max(D{j},1e-7);
    end
    res_d(:,i)=cell2mat(D);
    if i>1
        for j=1:length(D)
            D{j}=(D{j}*alpha+(1-alpha)*D_old{j});
        end  
    end
    
    figure(1)
    hold off
    plot(sqrt(sum((res_press(:,1:end-1)-res_press(:,2:end)).^2,1))./sqrt(sum((res_press(:,2:end)).^2,1)))
    hold on
    plot(sqrt(sum((res_d(:,1:end-1)-res_d(:,2:end)).^2,1))./sqrt(sum((res_d(:,2:end)).^2,1)))
    set(gca,'YScale','log')
    drawnow
    D_old=D;
    PRESSURE_old=PRESSURE;
    fprintf('alpha=%d\n',alpha);
alpha=min(alpha*1.01,0.8);
end