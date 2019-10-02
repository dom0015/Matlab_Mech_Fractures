function [ tmp1,tmp2,xx,yy,ftmp1,ftmp2,fxx,fyy,tmp1_node,tmp2_node] = streamlines_calc( u,node,elem ,gridx,gridy,fgridx,fgridy,fractures_positions)
%STREAMLINES_CALC Summary of this function goes here
%   Detailed explanation goes here
tmp_frac_dir=cell(length(fractures_positions),1);

for i=1:length(fractures_positions)
    dirx=fractures_positions{i}(1,1)-fractures_positions{i}(1,end);
    diry=fractures_positions{i}(2,1)-fractures_positions{i}(2,end);
    tmp_frac_dir{i}=[dirx*ones(length(fractures_positions{i}),1)  diry*ones(length(fractures_positions{i}),1)];
end

[frac_points,IA,IC]=unique(cell2mat(fractures_positions),'rows');
tmp=cell2mat(tmp_frac_dir);
frac_dir=1e1*tmp(IA,:);
nn=length(frac_points);

for i=1:3
x(:,i)=node(elem(:,i),1);
y(:,i)=node(elem(:,i),2);
z(:,i)=u(elem(:,i));
end
centers(:,1)=node(elem(:,1),1)+node(elem(:,2),1)+node(elem(:,3),1);
centers(:,2)=node(elem(:,1),2)+node(elem(:,2),2)+node(elem(:,3),2);
centers=centers/3;


grad_elem=[(z(:,1).*y(:,2) - z(:,2).*y(:,1) - z(:,1).*y(:,3) + z(:,3).*y(:,1) + z(:,2).*y(:,3) - z(:,3).*y(:,2))./(x(:,1).*y(:,2) - x(:,2).*y(:,1) - x(:,1).*y(:,3) + x(:,3).*y(:,1) + x(:,2).*y(:,3) - x(:,3).*y(:,2)), ...
                  -(z(:,1).*x(:,2) - z(:,2).*x(:,1) - z(:,1).*x(:,3) + z(:,3).*x(:,1) + z(:,2).*x(:,3) - z(:,3).*x(:,2))./(x(:,1).*y(:,2) - x(:,2).*y(:,1) - x(:,1).*y(:,3) + x(:,3).*y(:,1) + x(:,2).*y(:,3) - x(:,3).*y(:,2))];
            

x=gridx;
y=gridy;
[xx,yy]=meshgrid(x,y);
xx=xx(:);
yy=yy(:);
F_tmp=scatteredInterpolant(centers(:,1),centers(:,2),grad_elem(:,1),'linear','linear');
tmp1=F_tmp(xx,yy);
F_tmp=scatteredInterpolant(centers(:,1),centers(:,2),grad_elem(:,2),'linear','linear');
tmp2=F_tmp(xx,yy);


fx=fgridx;
fy=fgridy;
[fxx,fyy]=meshgrid(fx,fy);
%fxx=fxx(:);
%fyy=fyy(:);
F_tmp=scatteredInterpolant([centers(:,1);frac_points(:,1)],[centers(:,2);frac_points(:,2)],[grad_elem(:,1);frac_dir(:,1)],'linear','linear');
ftmp1=F_tmp(fxx,fyy);
F_tmp=scatteredInterpolant([centers(:,1);frac_points(:,1)],[centers(:,2);frac_points(:,2)],[grad_elem(:,2);frac_dir(:,2)],'linear','linear');
ftmp2=F_tmp(fxx,fyy);


A=sparse(elem(:),elem(:),reshape(repmat(grad_elem(:,1),3,1),length(elem)*3,1),length(node),length(node));
B=sparse(elem(:),elem(:),reshape(repmat(grad_elem(:,2),3,1),length(elem)*3,1),length(node),length(node));
C=sparse(elem(:),elem(:),ones(length(elem)*3,1),length(node),length(node));
tmp1_node=diag(A)./diag(C);
tmp2_node=diag(B)./diag(C);
end

