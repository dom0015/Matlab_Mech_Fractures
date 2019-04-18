%% 2d boundary value problem:
%   -div(k(x)*grad(u(x)))=f(x) inside of the domain
%   u(x)=U(x) at GammaD
%   k(x)*(du/dn)(x)=g(x) at GammaN
clear all
% close all

%% preparation of NODE, ELEM, bdFlag
height=1; width=1; h_elem=40; w_elem=40;
[ NODE,ELEM,bdFlag ] = rect_mesh( height, width, h_elem, w_elem );

%% add fractures
frac_start_end={[0.2 0.5], [0.8 0.5]
                [0.5 0.2], [0.5 0.8]};
[fractures, fractures_positions, no_fractures] = create_fractures( frac_start_end, NODE, h_elem );
[fractures_cell,fracture_matrice,intersections,lengths] = fracture2cells_geometry( fractures );
no_intersections = size(intersections,1);
[ NODE,ELEM,bdFlag,fractures_cell,fracture_matrice] = multi_fracture_tear( NODE,ELEM,fractures_cell,bdFlag,fracture_matrice);
%%
n_ELEM=length(ELEM);
n_NODE=length(NODE);
material_lambda=ones(n_ELEM,1);
material_mu=ones(n_ELEM,1);
F=[0*ones(n_NODE,1) 0*ones(n_NODE,1)];
MATERIAL=[material_lambda material_mu];

%%neumann 
% NBOUNDARY=[(1:(w_elem))' (2:(w_elem+1))'];
% right_side=[(n_NODE-(w_elem+1):n_NODE-1)' (n_NODE-(w_elem):n_NODE)'];
% NBOUNDARY=[NBOUNDARY; right_side];
% NVALUE=[0*ones(size(NBOUNDARY,1),1) 0*ones(size(NBOUNDARY,1),1)];
% NVALUE(1:end/2,:)=-NVALUE(1:end/2,:);
% 
% %% Neumann boundary at fracture sides
% % fracture 1
% NBOUNDARY=[NBOUNDARY; fracture_matrice{1}.above_nodes; fracture_matrice{1}.under_nodes];
% l=lengths(1)-1;
% temp=0.05*sin(linspace(0,5*pi,l))';
% NVALUE_fracture=[[0*ones(l,1) temp];[0*ones(l,1) 0*ones(l,1)]];
% NVALUE=[NVALUE; NVALUE_fracture];
% 
% % fracture 2
% NBOUNDARY=[NBOUNDARY; fracture_matrice{2}.above_nodes; fracture_matrice{2}.under_nodes];
% l=lengths(1)-1;
% temp=-0.05*sin(linspace(0,5*pi,l))';
% NVALUE_fracture=[[0*ones(l,1) temp];[0*ones(l,1) 0*ones(l,1)]];
% NVALUE=[NVALUE; NVALUE_fracture];

%% 1-bottom,2-right,3-top,4-left 
bdNeumann_x=double(bdFlag);
bdNeumann_x(bdFlag==1)=0;
bdNeumann_x(bdFlag==2)=0.2;
bdNeumann_x(bdFlag==3)=0;
bdNeumann_x(bdFlag==4)=0;
bdNeumann_y=double(0*bdFlag);
bdNeumann_y(bdFlag==3)=0.2;
bdNeumann_y(bdFlag==4)=0;
bdNeumann_val={bdNeumann_x,bdNeumann_y};

%% boundary conditions inputs
% temp=(NODE(:,2)==width);
% temp0=0*temp;
% side_up_x=[temp'; temp0']; side_up_x=side_up_x(:);
% side_up_y=[temp0'; temp']; side_up_y=side_up_y(:);
% temp=(NODE(:,2)==0);
% side_down_x=[temp'; temp0']; side_down_x=side_down_x(:);
% side_down_y=[temp0'; temp']; side_down_y=side_down_y(:);
% DBOUNDARY=side_down_x|side_down_y|side_up_x|side_up_y;
% temp=zeros(length(DBOUNDARY),1);
% temp(logical(side_up_x))=2.0;
% temp(logical(side_up_y))=1.0;
% DVALUE=temp(DBOUNDARY);
%% Dirichlet vlevo a dole
nahore=(NODE(:,2)==width)'; side_up_y=[nahore*0; nahore]; side_up_y=side_up_y(:);
vlevo=(NODE(:,1)==0)'; side_left_x=[vlevo; 0*vlevo]; side_left_x=side_left_x(:);
dole=(NODE(:,2)==0)'; side_down_y=[0*dole; dole]; side_down_y=side_down_y(:);
%vpravo=(NODE(:,1)==width)'; side_right_x=[vpravo; 0*vpravo]; side_right_x=side_right_x(:);

% right
DBOUNDARY=side_left_x|side_down_y;
temp=zeros(length(DBOUNDARY),1);
temp(logical(side_left_x))=0.0;
temp(logical(side_down_y))=0.0;
%temp(logical(side_right_x))=-0.1;

% up
% DBOUNDARY=side_left_x|side_down_y|side_up_y;
% temp=zeros(length(DBOUNDARY),1);
% temp(logical(side_left_x))=0.0;
% temp(logical(side_down_y))=0.0;
% temp(logical(side_up_y))=-0.1;

DVALUE=temp(DBOUNDARY);
FREENODE=true(n_NODE*2,1); FREENODE(DBOUNDARY)=false;

[A,b]=FEM_simple_fnc(NODE,ELEM,MATERIAL,F,bdFlag,bdNeumann_val);





%% modifications due to boundary conditions
u=zeros(n_NODE*2,1);
u(DBOUNDARY)=DVALUE;
b=b-A*u;
% for i=1:length(NVALUE)
%     x=NODE(NBOUNDARY(i,:),:);
%     len=norm(x(1,:)-x(2,:));
%     b(NBOUNDARY(i,:)*2-1)=b(NBOUNDARY(i,:)*2-1)+len*NVALUE(i,1)/2;
%     b(NBOUNDARY(i,:)*2)=b(NBOUNDARY(i,:)*2)+len*NVALUE(i,2)/2;
% end

%% solution of the resulting linear system and visualization
u(FREENODE)=A(FREENODE,FREENODE)\b(FREENODE);
coords1=NODE(:,1); coords2=NODE(:,2);
%hold on; triplot(ELEM,coords1,coords2,'g');
hold on; triplot(ELEM,coords1+u(1:2:end),coords2+u(2:2:end),'b');