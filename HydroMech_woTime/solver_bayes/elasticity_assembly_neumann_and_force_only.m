function [b]=elasticity_assembly_neumann_and_force_only(NODE,ELEM,bdFlagNeuman,bdNeuman_val,F)
%F=[0*ones(n_NODE,1) 0*ones(n_NODE,1)];
n_NODE=size(NODE,1);

coords1=NODE(:,1); coords2=NODE(:,2);
AREAS=polyarea(coords1(ELEM),coords2(ELEM),2);
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

bdFlag=bdFlagNeuman;
bdNeumann_x=bdNeuman_val{1};
bdNeumann_y=bdNeuman_val{2};
%% modifications due to boundary conditions
ELEM_Neuman_idx=find(sum(bdFlag,2)>0);
for i=ELEM_Neuman_idx'
    for j=1:3
        if(bdFlag(i,j)>0)
            el=ELEM(i,setdiff([1 2 3],mod(j+1,3)+1));
            x=NODE(el,:);
            len=norm(x(1,:)-x(2,:));
            val_x=bdNeumann_x(i,j);
            val_y=bdNeumann_y(i,j);
            b(el*2-1)=b(el*2-1)+len*val_x/2;
            b(el*2)=b(el*2)+len*val_y/2;
        end
    end
end