function [b]=elasticity_assembly_neumann_only(NODE,ELEM,bdFlagNeuman,bdNeuman_val)
%F=[0*ones(n_NODE,1) 0*ones(n_NODE,1)];
n_NODE=size(NODE,1);

b=zeros(2*n_NODE,1);

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