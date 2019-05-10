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
p=1e6; % pressure on Dir. b. c.
f = @(x)(0+0*x(:,1)+0*x(:,2)); % zatizeni
g_N=@(x)(0+0*x(:,1)+0*x(:,2)); % Neumann
mat_omega = 1e-15; % material - matrice
k = mat_omega*ones(length(elem),1);

d = exp([-6 -8.5]);
mat_frac = 1e-9*ones(no_fractures,1); % material - fractures
%alfa = 1e-7*ones(no_fractures,1); % prestup mezi trhlinami a matrici
alfa_inter = 1e-8*ones(no_intersections,1);

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
    Dirichlet_windows=[ 3   0.0     1.0     p
                        1   0.0     0.5     0
                        1   0.5     1.0     0];
    case 2
    Dirichlet_windows=[ 2   0.0     1.0     p
                        4   0.0     0.5     0
                        4   0.5     1.0     0];
    case 3
    Dirichlet_windows=[ 1   0.0     1.0     p
                        3   0.0     0.5     0
                        3   0.5     1.0     0];
    case 4
    Dirichlet_windows=[ 4   0.0     1.0     p
                        2   0.0     0.5     0
                        2   0.5     1.0     0];
    case 5
    Dirichlet_windows=[ 4   0    1    p
                        2   0     1     0];
end
Dirichlet_windows(:,2:3)=Dirichlet_windows(:,2:3)*10;

%% MATRICES ASSEMBLING
[u0, A, b, freeNode, downEdge, rightEdge, upEdge, leftEdge ] = FEM_windows( node, elem, h_elem, bdFlag, k, f, Dirichlet_windows, g_N );

[I,M]=interaction_matrix( fracture_matrice,node );
[F,b_f,u_f,freenode_f] = fractures_matrix( node,fracture_matrice,intersections,alfa_inter,lengths);
[B, freeNode, b, u0 ] = matrices_assembling( A, I, M, F, freeNode, freenode_f, b, b_f, u0, u_f );
u0(freeNode)=B(freeNode,freeNode)\b(freeNode);
%[u0(freeNode),flag,relres,iter,resvec]=gmres(B(freeNode,freeNode),b(freeNode),200,1e-6,1000);

%% extract flow
%% oblast slepena trhlinami

Q = extract_flow( A, u0(1:length(A)), node, elem, h_elem, Dirichlet_windows, downEdge, rightEdge, upEdge, leftEdge );
startx=10*ones(10,1);
starty=linspace(5,10,10);
[ tmp1,tmp2,xx,yy,ftmp1,ftmp2,fxx,fyy,node_fluxx,node_fluxy] = streamlines_calc( u0(1:length(A)),node,elem,linspace(0.1,9.9,23),linspace(0.1,9.9,23),...
    linspace(0,10,1000),linspace(0,10,1000),fractures_positions);

% Q = Q*mat_omega/h_elem;
disp(Q); 
%disp(sum(Q))

%% FIGURE
figure;
N=length(A);
h = trisurf(elem,node(:,1),node(:,2),u0(1:N));
alpha 0.8
h.EdgeColor = 'none';
colormap jet(1000)
colorbar
axis equal
hold on
view(0,90)
%set(gca,'ZDir','reverse')
quiver3(xx,yy,xx*0+p,-tmp1,-tmp2,tmp2*0,1.1,'Color','black','LineWidth',1.4)
%[x_g,y_g]=meshgrid(linspace(1,9,5),linspace(1,9,5));
n_sl=30;
x_g=[ones(n_sl,1);9*ones(n_sl,1)];
y_g=[linspace(0.1,9.9,n_sl)';linspace(0.1,9.9,n_sl)'];
frac_points=[x_g(:),y_g(:)];

hlines=streamline(fxx,fyy,ftmp1,ftmp2,frac_points(:,1),frac_points(:,2),[1e-2 1e7]);
set(hlines,'LineWidth',1.5,'Color','k');%[183 0 255]/256)

hlines=streamline(fxx,fyy,-ftmp1,-ftmp2,frac_points(:,1),frac_points(:,2),[1e-2 1e7]);
set(hlines,'LineWidth',1.5,'Color','k');%[183 0 255]/256)

%quiver(xx,yy,tmp1,tmp2,1.5,'Color','black','LineWidth',1.5)
%streamline(centers(:,1),centers(:,2),grad_elem(:,1),grad_elem(:,2),startx,starty)
%contour(node(:,1),node(:,2),u0(1:length(A)),10)
M = N;
for i=1:length(lengths)
    coord_x = fractures_positions{i}(:,1);
    coord_y = fractures_positions{i}(:,2);
    vals = u0(M+1:M+lengths(i));
    plot3(coord_x,coord_y,vals+1e6,'LineWidth',2); hold on
    diff_x = coord_x(2:end) - coord_x(1:end-1);
    diff_y = coord_y(2:end) - coord_y(1:end-1);
    steps_arg = sqrt(diff_x.^2 + diff_y.^2);
    steps_val = vals(2:end) - vals(1:end-1);
    %figure; plot(cumsum(steps_arg),steps_val)
    %figure; plot(cumsum(steps_arg),steps_arg)
    %figure; plot(cumsum(steps_arg),mat_fract(i)*steps_val./steps_arg)
    %figure; plot([0; cumsum(steps_arg)],vals)
    M=M+lengths(i);
end