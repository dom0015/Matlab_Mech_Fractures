function [Q,D,PRESSURE,ugrad,iter,x_elast,response_D] = solver_time_dependent(params,hydro_problem,problem_setting,SMALSE_params)
global const_cs_domain const_cs_fracture const_delta_t
hydro_problem.mat_frac=params;
no_fractures=hydro_problem.no_fractures;
lengths=hydro_problem.lengths;
d = 1e-4*ones(no_fractures,1);
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = d(i)*ones(lengths(i)-1,1);
end
eps_coupling=SMALSE_params.eps_coupling;
res_d=[];
res_press=[];
flag=0;
tic;
[PRESSURE__,u_old,ugrad,Q,PRESSURE]=tocouple_handle_modif(D,hydro_problem);
D_old=D;
PRESSURE_old=PRESSURE;
ugrad_old=ugrad;
elast_div=ugrad(:,1)*0; elast_div_diff=elast_div;
% ugrad_all=ugrad(:,1);
v_alpha=[];
kk=0;
rel_D_=[];
rel_P_=[];
response_D = D{1};
pressure_all = [];
u_old = 0*u_old;
load('u_old.mat')
for i=1:SMALSE_params.coupling_iter
    elast_div_old = elast_div;
    % one time step
    for j=1:100
        % hydro
        [PRESSURE__,u_new,ugrad,Q,PRESSURE,frac_grad,blocks]=tocouple_handle_modif_time...
            (D,hydro_problem,const_cs_domain,const_cs_fracture,u_old,elast_div_diff,const_delta_t);
        idx = [true(length(blocks.A),1); false(length(blocks.Au),1); false(length(blocks.F),1)];
        u_diff_norm = norm(u_old(idx)-u_new(idx));
        fprintf('TIME STEP i=%d: iteration=%d hydro diff=%d\n',i,j,u_diff_norm);
        u_old = u_new;
        pressure_all = [pressure_all PRESSURE__{1}];
        figure(13);
        N=length(hydro_problem.A);
        trisurf(hydro_problem.ELEMENTS,hydro_problem.POINTS(:,1),hydro_problem.POINTS(:,2),u_new(1:N),'LineStyle','none');
        view(0,90)
        colorbar
        colormap jet(1000)
        drawnow
%         pause()

        % mech
        [problem_setting] = feti.assembly_FETI_frac_rhs(problem_setting,PRESSURE,-ugrad);
        [D,problem_setting,x_elast] = smalse.SMALSE_solver(problem_setting,SMALSE_params);
        elast_div = elasticity_divergence(problem_setting,x_elast);
        elast_div_diff = elast_div - elast_div_old;
        temp = D{1}; temp(1)=temp(1)*2; temp(end)=temp(end)*2;
        response_D = [response_D temp];
        for k=1:length(D)
            D{k}=max(D{k},1e-10);  
        end
    end

    ugrad_all(:,i)=ugrad(:,1);
    tmp=cell2mat(PRESSURE__);
    res_press(:,i)=tmp(:);    
    res_d(:,i)=cell2mat(D);
    
    rel_P=sqrt(sum((res_press(:,1:end-1)-res_press(:,2:end)).^2,1))./sqrt(sum((res_press(:,2:end)).^2,1));
    rel_D=sqrt(sum((res_d(:,1:end-1)-res_d(:,2:end)).^2,1))./sqrt(sum((res_d(:,2:end)).^2,1));  
    
    % SMALSE precision to match coupled precision
    if ~isempty(rel_D)
        SMALSE_params.rel=1e-16;%min(rel_D(end),rel_P(end))/100000;
    end
    % coupling parameters update
    
    if flag==1
        kk=kk+1;
        rel_D_=[rel_D_ rel_D(end)];
        rel_P_=[rel_P_ rel_P(end)];
    end
  
    if SMALSE_params.print
        figure(12); 
        plot(linspace(3,7,length(u_new(N+1:end))),u_new(N+1:end))
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
        fprintf('i=%d: pres = %d || d = %d\n',i,rel_P(end),rel_D(end));
               
        if flag==1 && rel_P(end)<eps_coupling && rel_D(end)<eps_coupling
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

figure(111); imagesc(pressure_all); colorbar
end

