addpath(genpath('files_elasticity'))

%% BASIC GEOMETRY ---------------------------------------------------------
[POINTS,ELEMENTS,Dirichlet_boundaries,Neumann_boundaries,Neumann_normalx,Neumann_normaly,fractures,...
    fractures_positions,no_intersections,fractures_cell,fracture_matrice,fracture_elem_map] = ...
    create_geometry(Nxy,L1,L2,frac_start_end);

material_constants=100000*ones(length(ELEMENTS),2);
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
% N_bound(Neumann_boundaries==3)=1;
% N_bound(Neumann_boundaries==2)=1;

N_bound(Neumann_boundaries<0)=1;
%N_bound(Neumann_boundaries==-1.5)=1;

N_bound_value=cell(2,1);

tmp_nx=Neumann_normalx;
tmp_ny=Neumann_normaly;

for i=1:length(fracture_elem_map)
    [tmp_nx,tmp_ny] = fracture_boundary_values(fracture_elem_map{i},tmp_nx,tmp_ny,@(x)101*(1-4*(x-0.5).^2),@(x)101*(1-4*(x-0.5).^2));
end
N_bound_value{1}=tmp_nx;
N_bound_value{2}=tmp_ny;

%% FETI partitioning ------------------------------------------------------
[map,sub_nodes,sub_elem,sub_material_constants,sub_volume_force,sub_neumann,sub_neumann_val,sub_sizes,B_FETI,node_map_on,plot_func]=...
    partition_FETI(POINTS,ELEMENTS,material_constants,volume_force,N_bound,N_bound_value,sumbdomains_FETI);

%% Elasticity sub-Matrices assembly ---------------------------------------
A=cell(sumbdomains_FETI,1);
A_pinv=cell(sumbdomains_FETI,1);
A_null=cell(sumbdomains_FETI,1);
b=cell(sumbdomains_FETI,1);
for i=1:sumbdomains_FETI
    [ATemp,b{i}]=elasticity_assembly(sub_nodes{i},sub_elem{i},...
        sub_material_constants{i},sub_volume_force{i},...
        sub_neumann{i},sub_neumann_val{i});
    [APinvTemp,AKerTemp] = pinv_null(ATemp,100);
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
A=blkdiag(A{:});
b=cat(1,b{:});
A_plus = blkdiag(A_pinv{:});
R = blkdiag(A_null{:});
B_e=[B_dirichlet;B_FETI];
B_i = contact_inequalities(0*b,POINTS,fracture_matrice,node_map_on,size(A,1));
c_e=[B_dirichlet_rhs;zeros(size(B_FETI,1),1)];
c_i=zeros(size(B_i,1),1);

B_iupdate =@(x) contact_inequalities(x,POINTS,fracture_matrice,node_map_on,size(A,1));

rho=1e6;
n1=size(B_e,1);
n2=size(B_i,1);
O_e=sparse(n1,n1);
I_i=speye(n2,n2);
O_ie=sparse(n1,n2);

A_full=[A       B_e'    rho*B_i';...
        B_e    O_e     O_ie;
        rho*B_i O_ie'    -I_i];


b_full=[b;c_e;c_i];



%[x,flag,relres,iter,resvec] = gmres(A_full,b_full,100,1e-6,100);

tic;
x_full=A_full\b_full;
toc

B_i=B_iupdate(x_full);

%%
plot_func(x_full,fracture_matrice);
