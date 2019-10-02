function [g_P,M,NODE_above,NODE_under,crossed_L,crossed_R] = penalty_matrix_separated(u, POINTS, fracture_matrice,intersections)
%PENALTY_MATRIX Summary of this function goes here
%   Detailed explanation goes here
temp1=u(1:2:end);
temp2=u(2:2:end);
u_reshaped=[temp1 temp2]+POINTS;
u=u_reshaped'; u=u(:);

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

temp=(COORD_middles(2:end,:)-COORD_middles(1:end-1,:));
NORMALS_scaled=[temp(:,2) -temp(:,1)];
LENGTHS=sqrt(sum(temp.^2,2));
NORMALS=NORMALS_scaled./LENGTHS;
DISPLACEMENTS=COORD_above-COORD_under;

overlap_L=dot(NORMALS_scaled(1:end-1,:),DISPLACEMENTS,2);
overlap_R=dot(NORMALS_scaled(2:end,:),DISPLACEMENTS,2);
crossed_L=(overlap_L>0);
crossed_R=(overlap_R>0);

M=sparse(length(u),length(u));
g_P=sparse(length(u),1);

coords=[NORMALS_scaled(:,1) -NORMALS_scaled(:,1) NORMALS_scaled(:,2) -NORMALS_scaled(:,2)];
coords_L=coords(1:end-1,:);
coords_L(~crossed_L,:)=0;
coords_R=coords(2:end,:);
coords_R(~crossed_R,:)=0;
indices=[2*NODE_above-1 2*NODE_under-1 2*NODE_above 2*NODE_under];
for i=1:4 % rows
    for j=1:4 % columns
        temp_L=sparse(indices(:,i),indices(:,j),coords_L(:,i).*coords_L(:,j),length(u),length(u));
        temp_R=sparse(indices(:,i),indices(:,j),coords_R(:,i).*coords_R(:,j),length(u),length(u));
        M=M+temp_L+temp_R;
    end
    g_L=sparse(indices(:,i),1,coords_L(:,i).*overlap_L.*crossed_L,length(u),length(u));
    g_R=sparse(indices(:,i),1,coords_R(:,i).*overlap_R.*crossed_R,length(u),length(u));
    g_P=g_P+g_L+g_R;
end

M=M*2;
g_P=g_P*2;
end

