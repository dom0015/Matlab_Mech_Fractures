addpath(genpath('files_elasticity'))

%% PARAMETERS -------------------------------------------------------------
Nxy=201;
sumbdomains_FETI=ceil(Nxy/16)^2;
L1=1; L2=1;
frac_start_end={[0.2 0.4], [0.8 0.4]
    [0.2 0.7], [0.8 0.7]
     [0.2 0.2], [0.8 0.8]};
%     [0.3 0.2], [0.3 0.8]};

%% BASIC GEOMETRY ---------------------------------------------------------
nx=Nxy; ny=Nxy;
[coords1,coords2]=meshgrid(linspace(0,L1,nx),linspace(0,L2,ny));
coords1=coords1(:); coords2=coords2(:);
POINTS=[coords1,coords2];
ELEMENTS=rectangle_triangulation(nx,ny);
% NEUMANN - boundary edges
Neumann_boundaries=ELEMENTS*0;
Neumann_boundaries((1:((Nxy-1)*2):end),1)=1;
Neumann_boundaries((end-(Nxy-1)*2+1):2:end,2)=2;
Neumann_boundaries(((Nxy-1)*2):((Nxy-1)*2):end,2)=3;
Neumann_boundaries((2:2:(Nxy-1)*2),3)=4;
% DIRICHLET - boundary nodes
Dirichlet_boundaries=false(4,length(POINTS));
Dirichlet_boundaries(1,1:Nxy:end)=true;
Dirichlet_boundaries(2,(end-Nxy+1):end)=true;
Dirichlet_boundaries(3,Nxy:Nxy:end)=true;
Dirichlet_boundaries(4,1:Nxy)=true;
material_constants=ones(length(ELEMENTS),2);
volume_force=zeros(length(ELEMENTS),2);

%% Add fractures geometry -------------------------------------------------
[fractures, fractures_positions, no_fractures] = create_fractures( frac_start_end, POINTS, nx-1 );
[fractures_cell,fracture_matrice,intersections,lengths] = fracture2cells_geometry( fractures );
no_intersections = size(intersections,1);
[ POINTS,ELEMENTS,Neumann_boundaries,fractures_cell,fracture_matrice,Dirichlet_boundaries] = multi_fracture_tear( POINTS,ELEMENTS,fractures_cell ,Neumann_boundaries,fracture_matrice,Dirichlet_boundaries);

%% Boundary condition specification ---------------------------------------
u_0=zeros(2,length(POINTS));
D_bound=false(2,length(POINTS));
D_bound(2,Dirichlet_boundaries(1,:))=true;
D_bound(1,Dirichlet_boundaries(1,:))=true;
u_0(2,Dirichlet_boundaries(1,:))=0;
u_0(1,Dirichlet_boundaries(1,:))=0;
D_bound(2,Dirichlet_boundaries(3,:))=true;
u_0(2,Dirichlet_boundaries(3,:))=0.05;
D_bound=D_bound(:);
u_0=u_0(:);
N_bound=0*Neumann_boundaries;
N_bound_value=cell(2,1);
N_bound_value{1}=0*Neumann_boundaries;
N_bound_value{2}=0*Neumann_boundaries;

%% FETI partitioning ------------------------------------------------------
[map,sub_nodes,sub_elem,sub_material_constants,sub_volume_force,sub_neumann,sub_neumann_val,sub_sizes,B_FETI,node_map_on]=...
    partition_FETI(POINTS,ELEMENTS,material_constants,volume_force,N_bound,N_bound_value,sumbdomains_FETI);

%% Elasticity sub-Matrices assembly ---------------------------------------
A=cell(sumbdomains_FETI,1);
b=cell(sumbdomains_FETI,1);
for i=1:sumbdomains_FETI
    [A{i},b{i}]=elasticity_assembly(sub_nodes{i},sub_elem{i},...
                                    sub_material_constants{i},sub_volume_force{i},...
                                    sub_neumann{i},sub_neumann_val{i});
end

%% Matrix assembly --------------------------------------------------------
A_blokdiag=sparse(sub_sizes(end,3),sub_sizes(end,3));
b_blok=zeros(sub_sizes(end,3),1);
for i=1:sumbdomains_FETI
    A_blokdiag(sub_sizes(i,2):sub_sizes(i,3),sub_sizes(i,2):sub_sizes(i,3))=A{i};
    b_blok(sub_sizes(i,2):sub_sizes(i,3))=b{i};
end
dirichlet_nodes_old=find(D_bound);
dirichlet_nodes_new=node_map_on(dirichlet_nodes_old);
num_dnodes=length(dirichlet_nodes_new);
B_dirichlet=sparse(1:num_dnodes,dirichlet_nodes_new,ones(num_dnodes,1),num_dnodes,sub_sizes(end,3));
B_dirichlet_rhs=u_0(dirichlet_nodes_old);
B_full=[B_dirichlet;B_FETI];
B_rhs=[B_dirichlet_rhs;zeros(size(B_FETI,1),1)];
A_full=[A_blokdiag B_full';B_full sparse(size(B_full,1),size(B_full,1))];
b_full=[b_blok;B_rhs];
x_full=A_full\b_full;
x_elast=x_full(node_map_on);
x_elast1=POINTS(:,1)+x_elast(1:2:end);
x_elast2=POINTS(:,2)+x_elast(2:2:end);


%%
for i=1:size(frac_start_end,1)
    plot([frac_start_end{i,1}(1) frac_start_end{i,2}(1)],[frac_start_end{i,1}(2) frac_start_end{i,2}(2)],'k-','LineWidth',4)
end




figx=figure; triplot(ELEMENTS,POINTS(:,1),POINTS(:,2),'g');
hold on; triplot(ELEMENTS,x_elast1,x_elast2,'b');
hold off;
