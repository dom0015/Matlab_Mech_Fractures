function [u0,freeNode,B,b,elem,node,A,fractures_positions,lengths] = wrap_cg(mat_omega,mat_fracture,d,alfa_intersections,p_Dirichlet)


%% input parameters
% mat_omega = 1e-15; % material - matrice
% mat_fracture = 1e-9;
% d = exp([-6 -8.5]);
% alfa_intersections = 1e-8;
% p_Dirichlet=1e6; % pressure on Dir. b. c.


%% FEM with fractures
addpath(genpath('files'))
addpath(genpath('model_problems'))

%% GEOMETRY PARAMETERS
h_elem=40;
frac_start_end={[0.2 0.4], [0.8 0.4]
                [0.2 0.2], [0.8 0.8]};
            
%% GEOMETRY ASSEMBLING
[node,elem,bdFlag]=rect_mesh(10,10,h_elem,h_elem); % triangulace
[fractures, fractures_positions, no_fractures] = create_fractures( frac_start_end, node, h_elem );
[fractures_cell,fracture_matrice,intersections,lengths] = fracture2cells_geometry( fractures );
no_intersections = size(intersections,1);
[ node,elem ,bdFlag,fractures_cell,fracture_matrice] = multi_fracture_tear( node,elem,fractures_cell ,bdFlag,fracture_matrice);

%% MATERIAL AND BOUNDARY PARAMETERS
windows_type=5; % number of configuration of Dirichlet windows
f = @(x)(0+0*x(:,1)+0*x(:,2)); % zatizeni
g_N=@(x)(0+0*x(:,1)+0*x(:,2)); % Neumann
k = mat_omega*ones(length(elem),1);
mat_frac = mat_fracture*ones(no_fractures,1); % material - fractures
%alfa = 1e-7*ones(no_fractures,1); % prestup mezi trhlinami a matrici
alfa_inter = alfa_intersections*ones(no_intersections,1);

%% MATERIALS AND BOUNDARY ASSEMBLING
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = d(i)*ones(lengths(i)-1,1);
end

ALFA = cell(no_fractures,1);
MAT_FRAC = cell(no_fractures,1);
for i=1:no_fractures
    MAT_FRAC{i} = mat_frac(i).*D{i};
    ALFA{i} = mat_frac(i)./D{i};
end

[fracture_matrice] = fracture2cells_parameters( fracture_matrice,ALFA,MAT_FRAC );



%% Dirichletova okna ------------------------------------------------
% V kazdem radku:
% prvni hodnota - 1 dole; 2 vpravo; 3 nahore; 4 vlevo
% druha hodnota = a ... zacatek okna 
% treti hodnota = b ... konec okna ... 0 < a < b <= 1
% ctvrta hodnota - hodnota Dirichletovy podminky
switch windows_type
    case 1
    Dirichlet_windows=[ 3   0.0     1.0     p_Dirichlet
                        1   0.0     0.5     0
                        1   0.5     1.0     0];
    case 2
    Dirichlet_windows=[ 2   0.0     1.0     p_Dirichlet
                        4   0.0     0.5     0
                        4   0.5     1.0     0];
    case 3
    Dirichlet_windows=[ 1   0.0     1.0     p_Dirichlet
                        3   0.0     0.5     0
                        3   0.5     1.0     0];
    case 4
    Dirichlet_windows=[ 4   0.0     1.0     p_Dirichlet
                        2   0.0     0.5     0
                        2   0.5     1.0     0];
    case 5
    Dirichlet_windows=[ 4   0    1    p_Dirichlet
                        2   0     1     0];
end
Dirichlet_windows(:,2:3)=Dirichlet_windows(:,2:3)*10;

%% MATRICES ASSEMBLING
[u0, A, b, freeNode, downEdge, rightEdge, upEdge, leftEdge ] = FEM_windows( node, elem, h_elem, bdFlag, k, f, Dirichlet_windows, g_N );

[I,M]=interaction_matrix( fracture_matrice,node );
[F,b_f,u_f,freenode_f] = fractures_matrix( node,fracture_matrice,intersections,alfa_inter,lengths);
[B, freeNode, b, u0 ] = matrices_assembling( A, I, M, F, freeNode, freenode_f, b, b_f, u0, u_f );
% u0(freeNode)=B(freeNode,freeNode)\b(freeNode);