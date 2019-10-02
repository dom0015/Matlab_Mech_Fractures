addpath(genpath('files_hydro'))
%% MATERIAL AND BOUNDARY PARAMETERS
 % number of configuration of Dirichlet windows
p=1e6; % pressure on Dir. b. c.
f = @(x)(0+0*x(:,1)+0*x(:,2)); % zatizeni
g_N=@(x)(0+0*x(:,1)+0*x(:,2)); % Neumann


% %% GEOMETRY PARAMETERS - defined in elasticity
% L1=10; L2=10;
% Nxy=41;
% frac_start_end={[0.2 0.4], [0.8 0.4]
%                 [0.2 0.2], [0.8 0.8]
%                 [0.2 0.3], [0.8 0.3]};
%             
% %% GEOMETRY ASSEMBLING - defined in elasticity
% [POINTS,ELEMENTS,Dirichlet_boundaries,Neumann_boundaries,...
%     Neumann_normalx,Neumann_normaly,fractures,fractures_positions,...
%     no_intersections,fractures_cell,fracture_matrice,fracture_elem_map]...
%     = create_geometry(Nxy,L1,L2,frac_start_end);

%% ONLY FOR HYDRO
no_fractures=size(frac_start_end,1);
bdFlag=Neumann_boundaries;
bdFlag(bdFlag<0)=0;
bdFlag=bdFlag(:,[2,3,1]);
intersections=find_intersections(fractures);
lengths=zeros(1,no_fractures);
for i=1:no_fractures
    lengths(i)=length(fractures{i});
end
% [node_,elem_,bdFlag_]=rect_mesh(L1,L2,Nxy-1,Nxy-1); % triangulace
% [fractures_, fractures_positions_, no_fractures_] = create_fractures( frac_start_end, node, Nxy-1 );
% [fractures_cell_,fracture_matrice_,intersections_,lengths_] = fracture2cells_geometry( fractures );
% no_intersections = size(intersections,1);
%[ node,elem ,bdFlag,fractures_cell,fracture_matrice] = multi_fracture_tear( node,elem,fractures_cell ,bdFlag,fracture_matrice);
k = mat_omega_const*ones(length(ELEMENTS),1);
mat_frac = mat_frac_const*ones(no_fractures,1); % material - fractures
alfa_inter = alfa_inter_const*ones(no_intersections,1);

%% Dirichletova okna ------------------------------------------------
% V kazdem radku:
% prvni hodnota - 1 dole; 2 vpravo; 3 nahore; 4 vlevo
% druha hodnota = a ... zacatek okna 
% treti hodnota = b ... konec okna ... 0 < a < b <= 1
% ctvrta hodnota - hodnota Dirichletovy podminky
% switch windows_type
%     case 1
%     Dirichlet_windows=[ 3   0.0     1.0     p
%                         1   0.0     0.5     0
%                         1   0.5     1.0     0];
%     case 2
%     Dirichlet_windows=[ 2   0.0     1.0     p
%                         4   0.0     0.5     0
%                         4   0.5     1.0     0];
%     case 3
%     Dirichlet_windows=[ 1   0.0     1.0     p
%                         3   0.0     0.5     0
%                         3   0.5     1.0     0];
%     case 4
%     Dirichlet_windows=[ 4   0.0     1.0     p
%                         2   0.0     0.5     0
%                         2   0.5     1.0     0];
%     case 5
%     Dirichlet_windows=[ 2   0    1    p
%                         4   0    1     0];
% end
% Dirichlet_windows(:,2:3)=Dirichlet_windows(:,2:3);

Dirichlet_windows = D_windows( cislo_ulohy,p,n_windows,1/50 );

%% MATRICES ASSEMBLING
[u0, A, b, freeNode, downEdge, rightEdge, upEdge, leftEdge ] = FEM_windows( POINTS, ELEMENTS, Nxy-1, bdFlag, k, f, Dirichlet_windows, g_N );

hydro_problem.no_fractures=no_fractures;
hydro_problem.mat_frac=mat_frac;
hydro_problem.fracture_matrice=fracture_matrice;
hydro_problem.POINTS=POINTS;
hydro_problem.intersections=intersections;
hydro_problem.alfa_inter=alfa_inter;
hydro_problem.lengths=lengths;
hydro_problem.A=A;
hydro_problem.freeNode=freeNode;
hydro_problem.b=b;
hydro_problem.u0=u0;
hydro_problem.ELEMENTS=ELEMENTS;

hydro_problem.Dirichlet_windows=Dirichlet_windows;
hydro_problem.downEdge=downEdge;
hydro_problem.rightEdge=rightEdge;
hydro_problem.upEdge=upEdge;
hydro_problem.leftEdge=leftEdge;
hydro_problem.helem=Nxy-1;
