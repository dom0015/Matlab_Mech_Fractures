function [A,b]=elasticity_assembly(NODE,ELEM,MATERIAL,F,bdFlagNeuman,bdNeuman_val)

%MATERIAL 2 columns - lambda mu
%F - source - 2 columns
%%

material_lambda=MATERIAL(:,1);
material_mu=MATERIAL(:,2);
%F=[0*ones(n_NODE,1) 0*ones(n_NODE,1)];

coords1=NODE(:,1); coords2=NODE(:,2);
n_NODE=size(NODE,1);
n_ELEM=size(ELEM,1);

%% boundary conditions inputs
NBOUNDARY=[];
NVALUE=[];

%% other input data

c1111=material_lambda+2*material_mu;
c1122=material_lambda;
c1112=zeros(n_ELEM,1);
c2222=material_lambda+2*material_mu;
c2212=zeros(n_ELEM,1);
c1212=material_mu;
MATERIALS=[c1111 c1122 c1112 c1122 c2222 c2212 c1112 c2212 c1212];

%% construction of global "stiffness" matrix and rhs
AREAS=polyarea(coords1(ELEM),coords2(ELEM),2);
A=zeros(n_NODE*2);
% b=zeros(n_POINTS*2,1);
for i=1:n_ELEM
    % add local "stiffness" matrix
    x=NODE(ELEM(i,:),:);
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
    indices=reshape([ELEM(i,:)*2-1; ELEM(i,:)*2],1,6);
    A(indices,indices)=A(indices,indices)+A_local;
end

%% RHS vectorized
mid1x=(F(ELEM(:,2),1)+F(ELEM(:,3),1))/2; % midvalue
mid2x=(F(ELEM(:,1),1)+F(ELEM(:,3),1))/2; % midvalue
mid3x=(F(ELEM(:,1),1)+F(ELEM(:,2),1))/2; % midvalue
bt1x=AREAS.*(mid2x+mid3x)/6;
bt2x=AREAS.*(mid1x+mid3x)/6;
bt3x=AREAS.*(mid1x+mid2x)/6;
bx=accumarray(ELEM(:),[bt1x;bt2x;bt3x],[n_NODE 1]);

mid1y=(F(ELEM(:,2),2)+F(ELEM(:,3),2))/2; % midvalue
mid2y=(F(ELEM(:,1),2)+F(ELEM(:,3),2))/2; % midvalue
mid3y=(F(ELEM(:,1),2)+F(ELEM(:,2),2))/2; % midvalue
bt1y=AREAS.*(mid2y+mid3y)/6;
bt2y=AREAS.*(mid1y+mid3y)/6;
bt3y=AREAS.*(mid1y+mid2y)/6;
by=accumarray(ELEM(:),[bt1y;bt2y;bt3y],[n_NODE 1]);

b=zeros(2*n_NODE,1);
b(1:2:end)=bx;
b(2:2:end)=by;

%% modifications due to boundary conditions
% u=zeros(n_NODE*2,1);
% u(DBOUNDARY)=DVALUE;
% b=b-A*u;
% for i=1:length(NVALUE)
%     x=NODE(NBOUNDARY(i,:),:);
%     len=norm(x(1,:)-x(2,:));
%     b(NBOUNDARY(i,:)*2-1)=b(NBOUNDARY(i,:)*2-1)+len*NVALUE(i,1)/2;
%     b(NBOUNDARY(i,:)*2)=b(NBOUNDARY(i,:)*2)+len*NVALUE(i,2)/2;
% end
% 
% %% solution of the resulting linear system and visualization
% u(FREENODE)=A(FREENODE,FREENODE)\b(FREENODE);
% figure; triplot(ELEM,coords1,coords2,'g');
% hold on; triplot(ELEM,coords1+u(1:2:end),coords2+u(2:2:end),'b');