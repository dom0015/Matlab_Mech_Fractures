function [OVERLAP_L,OVERLAP_R] = calculate_overlap(u, POINTS, fracture_matrice,NODE_ABOVE,NODE_UNDER)
%CALCULATE_OVERLAP Summary of this function goes here
%   Detailed explanation goes here

temp1=u(1:2:end);
temp2=u(2:2:end);
u_reshaped=[temp1 temp2]+POINTS;

no_fractures=length(fracture_matrice);

OVERLAP_L=cell(no_fractures,1);
OVERLAP_R=cell(no_fractures,1);
%% normaly, smery otevreni, delky
for f=1:no_fractures
    NODE_above=NODE_ABOVE{f};
    NODE_under=NODE_UNDER{f};
    COORD_above=u_reshaped(NODE_above,:);
    COORD_under=u_reshaped(NODE_under,:);
    COORD_middles=(COORD_under+COORD_above)/2;

    temp=(COORD_middles(2:end,:)-COORD_middles(1:end-1,:));
    NORMALS_scaled=[temp(:,2) -temp(:,1)];
    DISPLACEMENTS=COORD_above-COORD_under;

    overlap_L=dot(NORMALS_scaled,DISPLACEMENTS(1:end-1,:),2);
    overlap_R=dot(NORMALS_scaled,DISPLACEMENTS(2:end,:),2);
    OVERLAP_L{f}=overlap_L;
    OVERLAP_R{f}=overlap_R;
end