function [u,A,b,FREENODE,fracture_matrice,POINTS,ELEMENTS,coords1,coords2,intersections] = elasticity(Nxy)
%ELASTICITY Summary of this function goes here
%   Detailed explanation goes here

%% example: preparation of a network
%  triangulation of a rectangular domain <0,L1>x<0,L2>
L1=1; L2=1;
nx=Nxy; ny=Nxy;
[coords1,coords2]=meshgrid(linspace(0,L1,nx),linspace(0,L2,ny));
coords1=coords1(:); coords2=coords2(:);
POINTS=[coords1,coords2];
ELEMENTS=rectangle_triangulation(nx,ny);

%% add fractures
frac_start_end={[0.2 0.4], [0.8 0.4]};
%                [0.2 0.2], [0.8 0.8]};
[fractures, fractures_positions, no_fractures] = create_fractures( frac_start_end, POINTS, nx-1 );
[fractures_cell,fracture_matrice,intersections,lengths] = fracture2cells_geometry( fractures );
no_intersections = size(intersections,1);
[ POINTS,ELEMENTS,bdFlag,fractures_cell,fracture_matrice] = multi_fracture_tear( POINTS,ELEMENTS,fractures_cell ,0,fracture_matrice);


coords1=POINTS(:,1); coords2=POINTS(:,2);
CENTERS_GRAVITY=[mean(coords1(ELEMENTS),2),mean(coords2(ELEMENTS),2)];
n_POINTS=size(POINTS,1);
n_ELEMENTS=size(ELEMENTS,1);

%% boundary conditions inputs
temp=(POINTS(:,2)==L2);
temp0=0*temp;
side_up_x=[temp'; temp0']; side_up_x=side_up_x(:);
side_up_y=[temp0'; temp']; side_up_y=side_up_y(:);
temp=(POINTS(:,2)==0);
side_down_x=[temp'; temp0']; side_down_x=side_down_x(:);
side_down_y=[temp0'; temp']; side_down_y=side_down_y(:);
DBOUNDARY=side_down_x|side_down_y|side_up_x|side_up_y;
temp=zeros(length(DBOUNDARY),1);
temp(logical(side_up_x))=0.0;
temp(logical(side_up_y))=0.0;
DVALUE=temp(DBOUNDARY);
FREENODE=true(n_POINTS*2,1); FREENODE(DBOUNDARY)=false;

NBOUNDARY=[(1:(ny-1))' (2:ny)'];
right_side=[(n_POINTS-ny:n_POINTS-1)' (n_POINTS-ny+1:n_POINTS)'];
NBOUNDARY=[NBOUNDARY; right_side];
NVALUE=[0*ones(size(NBOUNDARY,1),1) 0*ones(size(NBOUNDARY,1),1)];
NVALUE(1:floor(end/2),:)=-NVALUE(1:floor(end/2),:);

%% Neumann boundary at fracture sides
% fracture 1
NBOUNDARY=[NBOUNDARY; fracture_matrice{1}.above_nodes; fracture_matrice{1}.under_nodes];
l=lengths(1)-1;
% temp=0.05*sin(linspace(0,5*pi,l))';
temp=2*sin(linspace(0,5*pi,l))';
NVALUE_fracture=[[0*ones(l,1) temp];[0*ones(l,1) 0.5*ones(l,1)]];
NVALUE=[NVALUE; NVALUE_fracture];

% fracture 2
% NBOUNDARY=[NBOUNDARY; fracture_matrice{2}.above_nodes; fracture_matrice{2}.under_nodes];
% l=lengths(1)-1;
% temp=-0.05*sin(linspace(0,5*pi,l))';
% NVALUE_fracture=[[0*ones(l,1) temp];[0*ones(l,1) 0*ones(l,1)]];
% NVALUE=[NVALUE; NVALUE_fracture];


%% other input data
lambda=1;
mu=1;
c1111=lambda+2*mu;
c1122=lambda;
c1112=0;
c2222=lambda+2*mu;
c2212=0;
c1212=mu;
C=[c1111 c1122 c1112 c1122 c2222 c2212 c1112 c2212 c1212];
MATERIALS=repmat(C,n_ELEMENTS,1);
F=[0*ones(n_POINTS,1) 0*ones(n_POINTS,1)];

%% chosen area with different material
%temp=(CENTERS_GRAVITY(:,2)<0.5).*(CENTERS_GRAVITY(:,2)>0.2).*(CENTERS_GRAVITY(:,1)>0.5);
%MATERIALS(find(temp),:)=MATERIALS(find(temp),:)*4;

%% construction of global "stiffness" matrix and rhs
AREAS=polyarea(coords1(ELEMENTS),coords2(ELEMENTS),2);
%A=sparse(n_POINTS*2,n_POINTS*2);
if Nxy>121
    A=sparse([],[],[],n_POINTS*2,n_POINTS*2,n_POINTS*2*22);
else
    A=zeros(n_POINTS*2);
end
% b=zeros(n_POINTS*2,1);
for i=1:n_ELEMENTS
    % add local "stiffness" matrix
    x=POINTS(ELEMENTS(i,:),:);
    %     Bref=[-1.0  0.0 1.0 0.0 0.0 0.0
    %            0.0 -1.0 0.0 0.0 0.0 1.0
    %           -0.5 -0.5 0.0 0.5 0.5 0.0];
    Bref=[-1 1 0
        -1 0 1];
    DF=[x(2,1)-x(1,1) x(3,1)-x(1,1)
        x(2,2)-x(1,2) x(3,2)-x(1,2)];
    DFiT=[x(3,2)-x(1,2) x(1,2)-x(2,2)
        x(1,1)-x(3,1) x(2,1)-x(1,1)];
    DFiT=DFiT/det(DF);
    B=DFiT*Bref;
    BEeps=[B(1,1)  0  B(1,2)  0  B(1,3)  0
        0   B(2,1)  0  B(2,2)  0  B(2,3)
        reshape([B(2,:); B(1,:)],1,6)];
    CE=reshape(MATERIALS(i,:),3,3);
    A_local=BEeps'*CE*BEeps*AREAS(i);
    indices=reshape([ELEMENTS(i,:)*2-1; ELEMENTS(i,:)*2],1,6);
    A(indices,indices)=A(indices,indices)+A_local;
    % add local rhs
    %     FE=F(ELEMENTS(i,:),:)';
    %     b_local=FE(:)*AREAS(i);
    %     b(indices)=b(indices)+b_local;
end

%% RHS vectorized
mid1x=(F(ELEMENTS(:,2),1)+F(ELEMENTS(:,3),1))/2; % midvalue
mid2x=(F(ELEMENTS(:,1),1)+F(ELEMENTS(:,3),1))/2; % midvalue
mid3x=(F(ELEMENTS(:,1),1)+F(ELEMENTS(:,2),1))/2; % midvalue
bt1x=AREAS.*(mid2x+mid3x)/6;
bt2x=AREAS.*(mid1x+mid3x)/6;
bt3x=AREAS.*(mid1x+mid2x)/6;
bx=accumarray(ELEMENTS(:),[bt1x;bt2x;bt3x],[n_POINTS 1]);

mid1y=(F(ELEMENTS(:,2),2)+F(ELEMENTS(:,3),2))/2; % midvalue
mid2y=(F(ELEMENTS(:,1),2)+F(ELEMENTS(:,3),2))/2; % midvalue
mid3y=(F(ELEMENTS(:,1),2)+F(ELEMENTS(:,2),2))/2; % midvalue
bt1y=AREAS.*(mid2y+mid3y)/6;
bt2y=AREAS.*(mid1y+mid3y)/6;
bt3y=AREAS.*(mid1y+mid2y)/6;
by=accumarray(ELEMENTS(:),[bt1y;bt2y;bt3y],[n_POINTS 1]);

b=zeros(2*n_POINTS,1);
b(1:2:end)=bx;
b(2:2:end)=by;

%% modifications due to boundary conditions
u=zeros(n_POINTS*2,1);
u(DBOUNDARY)=DVALUE;
b=b-A*u;
for i=1:length(NVALUE)
    x=POINTS(NBOUNDARY(i,:),:);
    len=norm(x(1,:)-x(2,:));
    b(NBOUNDARY(i,:)*2-1)=b(NBOUNDARY(i,:)*2-1)+len*NVALUE(i,1)/2;
    b(NBOUNDARY(i,:)*2)=b(NBOUNDARY(i,:)*2)+len*NVALUE(i,2)/2;
end
A=sparse(A);
end

