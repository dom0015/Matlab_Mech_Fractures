function [ M1] = FEM1D_time(node,fracture)
%FEM1D Summary of this function goes here
%   Detailed explanation goes here
h1=sqrt(sum((node(fracture.above_nodes(:,1),:)-node(fracture.above_nodes(:,2),:)).^2,2));
h2=sqrt(sum((node(fracture.under_nodes(:,1),:)-node(fracture.under_nodes(:,2),:)).^2,2));
h=(h1+h2)/2;
n=length(h);

material=fracture.above_material*0+1;
tmp=material.*h;
tmp1=zeros(n+1,1);
tmp1(1:end-1)=tmp/6;
tmp2=zeros(n+1,1);
tmp2(1:end-1)=tmp/3;
tmp2(2:end)=tmp2(2:end)+tmp/3;
tmp3=zeros(n+1,1);
tmp3(2:end)=tmp/6;
M1=spdiags([tmp1 tmp2 tmp3],-1:1,n+1,n+1);

end

