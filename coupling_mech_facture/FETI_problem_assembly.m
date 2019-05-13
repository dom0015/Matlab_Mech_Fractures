addpath(genpath('files_elasticity'))
addpath(genpath('files_hydro'))
%% PARAMETERS -------------------------------------------------------------
Nxy=101;
L1=1; L2=1;
sumbdomains_FETI=ceil(Nxy/15)^2;

mat_const=10000*100000;
frac_press_val=1;
frac_start_end={[0.5 0.2], [0.9 0.2]
     [0.5 0.1], [0.9 0.5]
 %    [0.8 0.3], [0.8 0.9]
     [0.1 0.8], [0.9 0.8]
     [0.1 0.3], [0.7 0.9]
     [0.1 0.4], [0.7 0.4]
     [0.4 0.1], [0.4 0.5]
     [0.1 0.6], [0.7 0.6]};

%% BASIC GEOMETRY ---------------------------------------------------------
[POINTS,ELEMENTS,Dirichlet_boundaries,Neumann_boundaries,Neumann_normalx,Neumann_normaly,fractures,...
    fractures_positions,no_intersections,fractures_cell,fracture_matrice,fracture_elem_map] = ...
    create_geometry(Nxy,L1,L2,frac_start_end);

material_constants=mat_const*ones(length(ELEMENTS),2);
volume_force=zeros(length(ELEMENTS),2);
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
% N_bound(Neumann_boundaries==3)=1;
% N_bound(Neumann_boundaries==2)=1;

N_bound_fractures(Neumann_boundaries<0)=1;
N_bound_value{1}=Neumann_normalx;
N_bound_value{2}=Neumann_normaly;

%% FETI partitioning ------------------------------------------------------
[map,sub_nodes,sub_elem,sub_material_constants,sub_volume_force,...
    sub_neumann,sub_neumann_val,sub_sizes,B_FETI,node_map_on,plot_func,plot_func2,divide_neumann_boundary]=...
    partition_FETI(POINTS,ELEMENTS,material_constants,volume_force,N_bound,N_bound_value,sumbdomains_FETI);

%% Elasticity sub-Matrices assembly ---------------------------------------
A=cell(sumbdomains_FETI,1);
A_pinv=cell(sumbdomains_FETI,1);
A_null=cell(sumbdomains_FETI,1);
singular_values=zeros(sumbdomains_FETI,1);
b=cell(sumbdomains_FETI,1);
for i=1:sumbdomains_FETI
    [ATemp,b{i}]=elasticity_assembly(sub_nodes{i},sub_elem{i},...
        sub_material_constants{i},sub_volume_force{i},...
        sub_neumann{i},sub_neumann_val{i});
    [APinvTemp,AKerTemp,singular_values(i)] = pinv_null(ATemp,100);
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
b=cat(1,b{:});
A_plus = blkdiag(A_pinv{:});
R = blkdiag(A_null{:});
B_e=[B_dirichlet;B_FETI];
c_e=[B_dirichlet_rhs;zeros(size(B_FETI,1),1)];
c_i=0;
B_iupdate =@(x) contact_inequalities(x,POINTS,fracture_matrice,node_map_on,size(A,1));

%% create struct with problem parameters
problem_setting.A=A;
problem_setting.A_plus=A_plus;
problem_setting.B_e=B_e;
problem_setting.b=b;
problem_setting.R=R;
problem_setting.c_e=c_e;
problem_setting.c_i=c_i;
problem_setting.B_iupdate=B_iupdate;
problem_setting.fracture_matrice=fracture_matrice;
problem_setting.plot_func2=plot_func2;
problem_setting.mat_scale=mat_scale;

problem_setting.sub_nodes=sub_nodes;
problem_setting.sub_elem=sub_elem;
problem_setting.fracture_elem_map=fracture_elem_map;
problem_setting.divide_neumann_boundary=divide_neumann_boundary;
problem_setting.Neumann_normalx=Neumann_normalx;
problem_setting.Neumann_normaly=Neumann_normaly;
problem_setting.N_bound_fractures=N_bound_fractures;