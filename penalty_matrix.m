function [M,NODE_above,NODE_under,crossed,overlap] = penalty_matrix(u, POINTS, fracture_matrice,intersections)
%PENALTY_MATRIX Summary of this function goes here
%   Detailed explanation goes here
temp1=u(1:2:end);
temp2=u(2:2:end);
u_reshaped=[temp1 temp2]+POINTS;

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
COORD_middles=[u_reshaped(NODE_start,:); (COORD_under+COORD_above)/2; u_reshaped(NODE_end,:)];
%COORD_middles=POINTS([NODE_start; NODE_above; NODE_end],:);

M=sparse(length(u),length(u));
temp=(COORD_middles(3:end,:)-COORD_middles(1:end-2,:));

overlap=1/2*(temp(:,2).*(u_reshaped(NODE_above,1)-u_reshaped(NODE_under,1))...
            -temp(:,1).*(u_reshaped(NODE_above,2)-u_reshaped(NODE_under,2)));
crossed=(overlap>0);
%disp(overlap(crossed)')

coords=[temp(:,2) -temp(:,2) -temp(:,1) temp(:,1)];
coords(~crossed,:)=0;
indices=[2*NODE_above-1 2*NODE_under-1 2*NODE_above 2*NODE_under];
for i=1:4 % rows
    for j=1:4 % columns
        temp=sparse(indices(:,i),indices(:,j),coords(:,i).*coords(:,j),length(u),length(u));
        M=M+temp;
    end
end

end

