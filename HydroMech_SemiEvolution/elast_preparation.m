function [elast_problem,shared_data] = elast_preparation( elast_problem,shared_data )
frac_start_end = shared_data.frac_start_end;
Nxy = shared_data.Nxy;
L1 = shared_data.L1;
L2 = shared_data.L2;
par_tloustka_trhliny = elast_problem.par_tloustka_trhliny;
par_Lame_lambda = elast_problem.par_Lame_lambda;
par_Lame_mu = elast_problem.par_Lame_mu;
sumbdomains_FETI = elast_problem.sumbdomains_FETI;

no_fractures = size(frac_start_end,1);
fracture_direction=cell(1,no_fractures);
elast_problem.fracture_normal=cell(1,no_fractures);
for i=1:no_fractures
    if frac_start_end{i,1}(1) == frac_start_end{i,2}(1)
        fracture_direction{i} = 'v';
        elast_problem.fracture_normal{i} = [1; 0];
    elseif frac_start_end{i,1}(2) == frac_start_end{i,2}(2)
        fracture_direction{i} = 'h';
        elast_problem.fracture_normal{i} = [0; 1];
    else
        fracture_direction{i} = 'd';
        elast_problem.fracture_normal{i} = -[1; -1]/sqrt(2);
    end
end

d_trhlina=par_tloustka_trhliny*ones(no_fractures,1);

%% BASIC GEOMETRY ---------------------------------------------------------
[POINTS,ELEMENTS,Dirichlet_boundaries,Neumann_boundaries,Neumann_normalx,Neumann_normaly,fractures,...
    fractures_positions,no_intersections,fractures_cell,fracture_matrice,fracture_elem_map] = ...
    mesh.create_geometry(Nxy,L1,L2,frac_start_end);

for i=1:no_fractures
    temp = fracture_matrice{i}.above_nodes(:);
    idx_above = unique(temp(2:end-1));
    temp = fracture_matrice{i}.under_nodes(:);
    idx_under = unique(temp(2:end-1));
    if fracture_direction{i}=='h'
        POINTS(idx_above,2)=POINTS(idx_above,2)+d_trhlina(i)/2;
        POINTS(idx_under,2)=POINTS(idx_under,2)-d_trhlina(i)/2;
    elseif fracture_direction{i}=='v'
        POINTS(idx_above,1)=POINTS(idx_above,1)-d_trhlina(i)/2;
        POINTS(idx_under,1)=POINTS(idx_under,1)+d_trhlina(i)/2;
    else
        POINTS(idx_above,1)=POINTS(idx_above,1)-sqrt((d_trhlina(i)/2)^2/2);
        POINTS(idx_above,2)=POINTS(idx_above,2)+sqrt((d_trhlina(i)/2)^2/2);
        POINTS(idx_under,1)=POINTS(idx_under,1)+sqrt((d_trhlina(i)/2)^2/2);
        POINTS(idx_under,2)=POINTS(idx_under,2)-sqrt((d_trhlina(i)/2)^2/2);
    end
end

material_constants=[par_Lame_lambda*ones(length(ELEMENTS),1) par_Lame_mu*ones(length(ELEMENTS),1)];
volume_force=0*ones(length(ELEMENTS),2);

%% Boundary condition specification ---------------------------------------
u_0=zeros(2,length(POINTS));
D_bound=false(2,length(POINTS));
D_bound(2,Dirichlet_boundaries(1,:))=true;
D_bound(1,Dirichlet_boundaries(4,:))=true;
u_0(2,Dirichlet_boundaries(1,:))=0;
u_0(1,Dirichlet_boundaries(4,:))=0;
D_bound(1,Dirichlet_boundaries(2,:))=true;
u_0(1,Dirichlet_boundaries(2,:))=0;
D_bound(2,Dirichlet_boundaries(3,:))=true;
u_0(2,Dirichlet_boundaries(3,:))=0;
D_bound=D_bound(:);
u_0=u_0(:);
N_bound=0*Neumann_boundaries;
N_bound_fractures=0*Neumann_boundaries;
N_bound_fractures(Neumann_boundaries<0)=1;
N_bound_value{1}=Neumann_normalx;
N_bound_value{2}=Neumann_normaly;

%% FETI partitioning ------------------------------------------------------
[map,sub_nodes,sub_elem,sub_material_constants,sub_volume_force,...
    sub_neumann,sub_neumann_val,sub_sizes,B_FETI,node_map_on,plot_func,plot_func2,divide_neumann_boundary,divide_neumann_boundary_and_force]=...
    feti.partition_FETI(POINTS,ELEMENTS,material_constants,volume_force,N_bound,N_bound_value,sumbdomains_FETI);

%% Elasticity sub-Matrices assembly ---------------------------------------
A=cell(sumbdomains_FETI,1);
A_pinv=cell(sumbdomains_FETI,1);
A_null=cell(sumbdomains_FETI,1);
singular_values=zeros(sumbdomains_FETI,1);
b=cell(sumbdomains_FETI,1);
NAPETI=cell(sumbdomains_FETI,1);
for i=1:sumbdomains_FETI
    [ATemp,b{i},NAPETI{i}]=a_ela.elasticity_assembly(sub_nodes{i},sub_elem{i},...
        sub_material_constants{i},sub_volume_force{i},...
        sub_neumann{i},sub_neumann_val{i});
    [APinvTemp,AKerTemp,singular_values(i)] = feti.pinv_null(ATemp,100000);
    A_pinv{i}=sparse(APinvTemp);
    A_null{i}=sparse(AKerTemp);
    A{i}=sparse(ATemp);
end

%% Dirichlet boundary assembly
dirichlet_nodes_old=find(D_bound);
dirichlet_nodes_new=node_map_on(dirichlet_nodes_old);
num_dnodes=length(dirichlet_nodes_new);
B_dirichlet=sparse(1:num_dnodes,dirichlet_nodes_new,ones(num_dnodes,1),num_dnodes,sub_sizes(end,3));
B_dirichlet_rhs=u_0(dirichlet_nodes_old);

%% Matrix assembly --------------------------------------------------------
mat_scale=max(singular_values);
A=blkdiag(A{:});
NAPETI=blkdiag(NAPETI{:});
b=cat(1,b{:});
A_plus = blkdiag(A_pinv{:});
R = blkdiag(A_null{:});
B_e=[B_dirichlet;B_FETI];
c_e=[B_dirichlet_rhs;zeros(size(B_FETI,1),1)];
c_i=0;
B_iupdate =@(x) a_ela.contact_inequalities(x,POINTS,fracture_matrice,node_map_on,size(A,1));

%% create struct with problem parameters
elast_problem.A=A;
elast_problem.A_plus=A_plus;
elast_problem.B_e=B_e;
elast_problem.b=b;
elast_problem.R=R;
elast_problem.c_e=c_e;
elast_problem.c_i=[];
elast_problem.B_iupdate=B_iupdate;
elast_problem.plot_func2=plot_func2;
elast_problem.mat_scale=mat_scale;

elast_problem.sub_nodes=sub_nodes;
elast_problem.sub_elem=sub_elem;
elast_problem.fracture_elem_map=fracture_elem_map;
elast_problem.divide_neumann_boundary=divide_neumann_boundary;
elast_problem.divide_neumann_boundary_and_force=divide_neumann_boundary_and_force;
elast_problem.Neumann_normalx=Neumann_normalx;
elast_problem.Neumann_normaly=Neumann_normaly;
elast_problem.N_bound_fractures=N_bound_fractures;
elast_problem.NAPETI=NAPETI;
elast_problem.map=map;
elast_problem.node_map_on=node_map_on;
elast_problem.fracture_matrice=fracture_matrice;
elast_problem.B_i=[];
elast_problem.lambda_ker=[];

shared_data.Neumann_boundaries=Neumann_boundaries;
shared_data.fractures=fractures;
shared_data.POINTS=POINTS;
shared_data.ELEMENTS=ELEMENTS;
shared_data.no_intersections=no_intersections;
shared_data.fracture_matrice=fracture_matrice;