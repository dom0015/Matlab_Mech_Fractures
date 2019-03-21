%% 2d boundary value problem:
%   -div(k(x)*grad(u(x)))=f(x) inside of the domain
%   u(x)=U(x) at GammaD
%   k(x)*(du/dn)(x)=g(x) at GammaN

%% example: preparation of a network
%  triangulation of a rectangular domain <0,L1>x<0,L2>
L1=1; L2=1;
nx=41; ny=41;
[coords1,coords2]=meshgrid(linspace(0,L1,nx),linspace(0,L2,ny));
coords1=coords1(:); coords2=coords2(:);
POINTS=[coords1,coords2];
ELEMENTS=rectangle_triangulation(nx,ny);

%% add fractures
frac_start_end={[0.2 0.4], [0.8 0.4]
                [0.2 0.2], [0.8 0.8]};
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
NVALUE(1:end/2,:)=-NVALUE(1:end/2,:);

%% Neumann boundary at fracture sides
% fracture 1
NBOUNDARY=[NBOUNDARY; fracture_matrice{1}.above_nodes; fracture_matrice{1}.under_nodes];
l=lengths(1)-1;
temp=0.05*sin(linspace(0,5*pi,l))';
NVALUE_fracture=[[0*ones(l,1) temp];[0*ones(l,1) 0*ones(l,1)]];
NVALUE=[NVALUE; NVALUE_fracture];

% fracture 2
NBOUNDARY=[NBOUNDARY; fracture_matrice{2}.above_nodes; fracture_matrice{2}.under_nodes];
l=lengths(1)-1;
temp=-0.05*sin(linspace(0,5*pi,l))';
NVALUE_fracture=[[0*ones(l,1) temp];[0*ones(l,1) 0*ones(l,1)]];
NVALUE=[NVALUE; NVALUE_fracture];

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
A=zeros(n_POINTS*2);
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

%% solution of the resulting linear system and visualization
u(FREENODE)=A(FREENODE,FREENODE)\b(FREENODE);
figure; triplot(ELEMENTS,coords1,coords2,'g');
hold on; triplot(ELEMENTS,coords1+u(1:2:end),coords2+u(2:2:end),'b');

temp1=u(1:2:end);
temp2=u(2:2:end);
u_reshaped=[temp1 temp2]+POINTS;

%% normaly, smery otevreni, delky
% ELEM_above=fracture_matrice{1}.above_nodes;
% orig1=POINTS(ELEM_above(:,1),:);
% orig2=POINTS(ELEM_above(:,2),:);
% temp=orig2-orig1;
% LEN_frac=sqrt(sum(temp.^2,2));                          % delky
% NORM_frac=[temp(:,2) -temp(:,1)]./[LEN_frac LEN_frac];   % normaly
% above1=u_reshaped(ELEM_above(:,1),:)-orig1;
% above2=u_reshaped(ELEM_above(:,2),:)-orig2;
% ELEM_under=fracture_matrice{1}.under_nodes;
% under1=u_reshaped(ELEM_under(:,1),:)-orig1;
% under2=u_reshaped(ELEM_under(:,2),:)-orig2;
% overlap=dot(NORM_frac,above1,2)+dot(NORM_frac,above2,2)...
%     -dot(NORM_frac,under1,2)-dot(NORM_frac,under2,2);
% overlap=overlap/2;
% crossed=(overlap<0);
% d_overlap=0*u_reshaped;
% d_overlap(ELEM_above(crossed,1),:)=d_overlap(ELEM_above(crossed,1),:)+NORM_frac(crossed,:);
% d_overlap(ELEM_above(crossed,2),:)=d_overlap(ELEM_above(crossed,2),:)+NORM_frac(crossed,:);
% d_overlap(ELEM_under(crossed,1),:)=d_overlap(ELEM_under(crossed,1),:)-NORM_frac(crossed,:);
% d_overlap(ELEM_under(crossed,2),:)=d_overlap(ELEM_under(crossed,2),:)-NORM_frac(crossed,:);
% d_overlap=d_overlap/2;
% penalty_sqrt=-sum(overlap(crossed));
% penalty=penalty_sqrt^2;
% d_penalty=2*penalty_sqrt*(-d_overlap);

%% normaly, smery otevreni, delky
ELEM_above=fracture_matrice{1}.above_nodes;
ELEM_under=fracture_matrice{1}.under_nodes;
NODE_start=ELEM_above(1,1);
NODE_end=ELEM_above(end,2);
NODE_above=ELEM_above(2:end,1);
NODE_under=ELEM_under(2:end,1);
if ~isempty(intersections)
    NODE_above_2=ELEM_above(1:end-1,2);
    temp=(NODE_above~=NODE_above_2); % temp is nonzero in case of crossing fractures
    shift=0;
    for i=find(temp)'
        NODE_above=[NODE_above(1:(i+shift)); NODE_above_2(i); NODE_above((i+shift+1):end)];
        NODE_under=[NODE_under(1:(i+shift)); NODE_above_2(i); NODE_under((i+shift+1):end)];
        shift=shift+1;
    end
end
COORD_above=u_reshaped(NODE_above,:);
COORD_under=u_reshaped(NODE_under,:);
%COORD_middles=[u_reshaped(NODE_start,:); (COORD_under+COORD_above)/2; u_reshaped(NODE_end,:)];
COORD_middles=POINTS([NODE_start; NODE_above; NODE_end],:);
temp=COORD_middles(1:end-1,:)-COORD_middles(2:end,:);
temp=[temp(:,2) -temp(:,1)];
NORM_frac=(temp(1:end-1,:)+temp(2:end,:))/2; % averaged normal vector (for inside nodes)
%COORD_middles=COORD_middles(2:end-1,:);
VECT_above=COORD_above-COORD_middles(2:end-1,:);
VECT_under=COORD_under-COORD_middles(2:end-1,:);
overlap=dot(NORM_frac,VECT_above,2)-dot(NORM_frac,VECT_under,2);
crossed=(overlap<0);
% i_xa=repmat(2*NODE_above-1,4,1);
% i_ya=repmat(2*NODE_above,4,1);
% i_xu=repmat(2*NODE_under-1,4,1);
% i_yu=repmat(2*NODE_under,4,1);
% i=[i_xa; i_ya; i_xu; i_yu];
% j=[2*NODE_above-1; 2*NODE_above; 2*NODE_under-1; 2*NODE_under];
% j=repmat(j,4,1);
M=sparse(length(u),length(u));
temp=(COORD_middles(3:end,:)-COORD_middles(1:end-2,:));
coord_A=temp(:,1);
coord_B=temp(:,2);
coords=[temp(:,2) -temp(:,2) -temp(:,1) temp(:,1)];
coords(~crossed,:)=0;
indices=[2*NODE_above-1 2*NODE_under-1 2*NODE_above 2*NODE_under];
for i=1:4 % rows
    for j=1:4 % columns
        temp=sparse(indices(:,i),indices(:,j),coords(:,i).*coords(:,j),length(u),length(u));
        M=M+temp;
    end
end

plot(u_reshaped(NODE_above(crossed),1),u_reshaped(NODE_above(crossed),2),'r*')
plot(u_reshaped(NODE_under(crossed),1),u_reshaped(NODE_under(crossed),2),'m*')
figure; plot(overlap);

%% find overlapping nodes
% for p=1:length(u_reshaped)
%     s=u_reshaped(p,:);
%     for e=1:length(ELEMENTS)
%         t=ELEMENTS(e,:);
%         if p==t(1) || p==t(2) || p==t(3)
%             continue
%         end
%         a=u_reshaped(t(1),:);
%         b=u_reshaped(t(2),:);
%         c=u_reshaped(t(3),:);
%         if point_inside_triangle(s,a,b,c)
%             plot(s(1),s(2),'*r')
%         end
%     end
% end

