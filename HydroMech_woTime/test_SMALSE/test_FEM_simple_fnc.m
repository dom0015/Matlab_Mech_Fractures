%% 2d boundary value problem:
%   -div(k(x)*grad(u(x)))=f(x) inside of the domain
%   u(x)=U(x) at GammaD
%   k(x)*(du/dn)(x)=g(x) at GammaN
clear all
close all

%% preparation of NODE, ELEM, bdFlag
height=1; width=1; h_elem=40; w_elem=40;
[ NODE,ELEM,bdFlag ] = rect_mesh( height, width, h_elem, w_elem );

%% add fractures
frac_start_end={[0.2 0.4], [0.8 0.4]
                [0.2 0.2], [0.8 0.8]};
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

[A,b]=FEM_simple_fnc(NODE,ELEM,MATERIAL,F,[],[]);


%% boundary conditions inputs
temp=(NODE(:,2)==width);
temp0=0*temp;
side_up_x=[temp'; temp0']; side_up_x=side_up_x(:);
side_up_y=[temp0'; temp']; side_up_y=side_up_y(:);
temp=(NODE(:,2)==0);
side_down_x=[temp'; temp0']; side_down_x=side_down_x(:);
side_down_y=[temp0'; temp']; side_down_y=side_down_y(:);
DBOUNDARY=side_down_x|side_down_y|side_up_x|side_up_y;
temp=zeros(length(DBOUNDARY),1);
temp(logical(side_up_x))=0.0;
temp(logical(side_up_y))=0.0;
DVALUE=temp(DBOUNDARY);
FREENODE=true(n_NODE*2,1); FREENODE(DBOUNDARY)=false;

NBOUNDARY=[];
NVALUE=[];
%% modifications due to boundary conditions
u=zeros(n_NODE*2,1);
u(DBOUNDARY)=DVALUE;
b=b-A*u;
for i=1:length(NVALUE)
    x=NODE(NBOUNDARY(i,:),:);
    len=norm(x(1,:)-x(2,:));
    b(NBOUNDARY(i,:)*2-1)=b(NBOUNDARY(i,:)*2-1)+len*NVALUE(i,1)/2;
    b(NBOUNDARY(i,:)*2)=b(NBOUNDARY(i,:)*2)+len*NVALUE(i,2)/2;
end

%% solution of the resulting linear system and visualization
u(FREENODE)=A(FREENODE,FREENODE)\b(FREENODE);
coords1=NODE(:,1); coords2=NODE(:,2);
figure; triplot(ELEM,coords1,coords2,'g');
hold on; triplot(ELEM,coords1+u(1:2:end),coords2+u(2:2:end),'b');