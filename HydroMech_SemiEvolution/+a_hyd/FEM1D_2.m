function [ A,b,u,freenode ] = FEM1D_2(node,fracture)
%FEM1D Summary of this function goes here
%   Detailed explanation goes here
h1=sqrt(sum((node(fracture.above_nodes(:,1),:)-node(fracture.above_nodes(:,2),:)).^2,2));
h2=sqrt(sum((node(fracture.under_nodes(:,1),:)-node(fracture.under_nodes(:,2),:)).^2,2));
h=(h1+h2)/2;
material=fracture.material;
n=length(h);
tmp=material./h;

tmp1=zeros(n+1,1);
tmp1(1:end-1)=-tmp;
tmp2=zeros(n+1,1);
tmp2(1:end-1)=tmp;
tmp2(2:end)=tmp2(2:end)+tmp;
tmp3=zeros(n+1,1);
tmp3(2:end)=-tmp;
A=spdiags([tmp1 tmp2 tmp3],-1:1,n+1,n+1);
freenode=true(n+1,1);
u=zeros(n+1,1);
if fracture.left_boundary(1)==1
    freenode(1)=0;
    u(1)=left_boundary(2);
end

if fracture.right_boundary(1)==1
    freenode(end)=0;
    u(end)=right_boundary(2);
end

b=-A*u;
end

