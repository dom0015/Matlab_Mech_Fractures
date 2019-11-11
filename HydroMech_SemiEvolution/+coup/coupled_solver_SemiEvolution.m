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
res_d=[];
res_press=[];
alpha=1;
flag=0;
tic;
[PRESSURE__,u_old,ugrad,Q,PRESSURE]=coup.hydro_Stationary(D,hydro_problem);
% elast_div=ugrad(:,1)*0; elast_div_diff=elast_div;
kk=0;
rel_D_=[];
rel_P_=[];
kkk=[];
response_D = cell2mat(D);
%u_old = 0*u_old;
for i=1:SMALSE_params.coupling_iter
    ugrad_old=ugrad;
    [PRESSURE__,u0_,ugrad,Q,PRESSURE,frac_grad]=coup.hydro_basicSemiEvolution...
        (D,hydro_problem,u_old);%,elast_div_diff);
    
    if i>1
        max(abs(ugrad_old-ugrad))
        kkk=[kkk max(abs(ugrad_old-ugrad))];
    end
    u_old = u0_;
    tmp=cell2mat(PRESSURE__);
    res_press(:,i)=tmp(:);
    
    [elast_problem] = feti.assembly_FETI_frac_rhs(elast_problem,PRESSURE,-ugrad*par_BiotWillis);
    [D,elast_problem,x_elast] = smalse.SMALSE_solver(elast_problem,SMALSE_params);
    %     elast_div_old = elast_div;
    %     elast_div = elasticity_divergence(elast_problem,x_elast,SMALSE_params.print);
    %     elast_div_diff = (elast_div - elast_div_old)*par_BiotWillis;
    response_D = [response_D cell2mat(D)];
    for j=1:length(D)
        D{j}=max(D{j},0)+par_a0;
    end
    
    res_d(:,i)=cell2mat(D);
    
    if i>1
        fprintf('----%d----',max((max(res_d(:,i)./res_d(:,i-1),res_d(:,i-1)./res_d(:,i))))-1)
    end
    rel_P=sqrt(sum((res_press(:,1:end-1)-res_press(:,2:end)).^2,1))./sqrt(sum((res_press(:,2:end)).^2,1));
    rel_D=sqrt(sum((res_d(:,1:end-1)-res_d(:,2:end)).^2,1))./sqrt(sum((res_d(:,2:end)).^2,1));
    
    % SMALSE precision to match coupled precision
    %     if ~isempty(rel_D)
    %         SMALSE_params.rel=1e-16;%min(rel_D(end),rel_P(end))/100000;
    %     end
    
    if flag==1
        kk=kk+1;
        rel_D_=[rel_D_ rel_D(end)];
        rel_P_=[rel_P_ rel_P(end)];
    end
    
    if SMALSE_params.print
        figure(105);
        N=length(hydro_problem.A);
        trisurf(hydro_problem.ELEMENTS,hydro_problem.POINTS(:,1),hydro_problem.POINTS(:,2),u0_(1:N),'LineStyle','none');
        colormap jet; axis equal; view(0,90); title('Hydraulic pressure')
        figure(106);
        %         plot(linspace(3,7,length(u0_(N+1:end))),u0_(N+1:end))
        plot(u0_(N+1:end))
        grid on
        title("Hydro N+1:end");
        figure(107);
        %         tmp=frac_grad{1};
        %         tmp(1)=tmp(1)/2;
        %         tmp(end)=tmp(end)/2;
        %         plot(linspace(3,7,length(-1e-5*tmp)),-1e-5*tmp)
        tmp=cell2mat(frac_grad);
        plot(tmp)
        title("Fracture gradient")
    end
    
    if ~isempty(rel_D)
        fprintf('i=%d: alpha = %.4f , beta = %.2f || pres = %d || d = %d\n',i,alpha,flag,rel_P(end),rel_D(end));
        
        if flag==1 && alpha>=0.00001 && rel_P(end)<eps_coupling && rel_D(end)<eps_coupling
            break
        end
        if alpha <0.00001
            break
        end
    end
    
    %     D_old=D;
    %     PRESSURE_old=PRESSURE;
    %     ugrad_old=ugrad;
    %
    flag=1;
end

if SMALSE_params.print_couple
    figure(201)
    subplot(1,2,1)
    hold on
    plot(rel_P_,'k-','LineWidth',1)
    set(gca,'YScale','log')
    title("rel P");
    figure(201)
    subplot(1,2,2)
    yyaxis right
    hold on
    plot(rel_D_,'k-','LineWidth',1)
    set(gca,'YScale','log')
    title("rel D");
end
iter=length(rel_D_);
end

