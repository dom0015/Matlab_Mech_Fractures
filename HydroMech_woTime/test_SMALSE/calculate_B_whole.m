function [B_all] = calculate_B_whole(u, NODE, fracture_matrice,intersections)
%PENALTY_MATRIX Summary of this function goes here
%   Detailed explanation goes here

temp1=u(1:2:end);
temp2=u(2:2:end);
u_reshaped=[temp1 temp2]+NODE;
u=u_reshaped'; u=u(:);

%% local matrices
Z=zeros(4);
% M_R_variable=[0 0 -0.5 0.5; 0 0 -0.5 0.5; 0   0    0.5 -0.5;  0    0   0.5 -0.5];
M_R_constant=[0 0  0   0  ; 0 0  0   0  ; 0.5 0.5 -0.5 -0.5; -0.5 -0.5 0.5  0.5];
% M_L_variable=[-0.5 0.5  0    0  ; -0.5  0.5 0   0  ; 0.5 -0.5 0 0; 0.5 -0.5 0 0];
M_L_constant=[ 0.5 0.5 -0.5 -0.5; -0.5 -0.5 0.5 0.5; 0    0   0 0; 0    0   0 0];
M_R=M_R_constant;%+M_R_variable;
M_L=M_L_constant;%+M_L_variable;
M_R=[Z -M_R; M_R Z];
M_L=[Z -M_L; M_L Z];

no_fractures=length(fracture_matrice);
B=cell(no_fractures,1);
%% normaly, smery otevreni, delky
for f=1:no_fractures
    ELEM_above=fracture_matrice{f}.above_nodes;
    ELEM_under=fracture_matrice{f}.under_nodes;
    NODE_start=ELEM_above(1,1);
    NODE_end=ELEM_above(end,2);
    NODE_above=ELEM_above(2:end,1);
    NODE_under=ELEM_under(2:end,1);
    if ~isempty(intersections)
        NODE_above_2=ELEM_above(1:end-1,2);
        temp_a=(NODE_above~=NODE_above_2); % temp is nonzero in case of intersecting fractures
        NODE_under_2=ELEM_under(1:end-1,2);
%         temp_u=(NODE_under~=NODE_under_2); % temp is nonzero in case of intersecting fractures
        shift=0;
        for i=find(temp_a)'
            NODE_above=[NODE_above(1:(i+shift-1)); NODE_above_2(i); NODE_above((i+shift):end)];
            NODE_under=[NODE_under(1:(i+shift-1)); NODE_under_2(i); NODE_under((i+shift):end)];
            shift=shift+1;
        end
    end
    NODE_above = [NODE_start; NODE_above; NODE_end];
    NODE_under = [NODE_start; NODE_under; NODE_end];

    ua_L=NODE_above(1:end-1,:);
    uu_L=NODE_under(1:end-1,:);
    ua_R=NODE_above(2:end,:);
    uu_R=NODE_under(2:end,:);
    temp=[ua_L uu_L ua_R uu_R];
    indices=[2*temp-1 2*temp];

    N=length(ua_L)+length(uu_L);
    mat_B=sparse(N,length(u));
    i_B=1;

    for i=1:length(ua_L)
        idx=indices(i,:);
        Hd_L_u=M_L*u(idx);
        Hd_R_u=M_R*u(idx);
%         mat_B(i_B,idx)=Hd_L_u; i_B=i_B+1;
%         mat_B(i_B,idx)=Hd_R_u; i_B=i_B+1;
        mat_B=mat_B+sparse(i_B,idx,Hd_L_u,N,length(u)); i_B=i_B+1;
        mat_B=mat_B+sparse(i_B,idx,Hd_R_u,N,length(u)); i_B=i_B+1;
    end
    B{f}=mat_B;
end
B_all=B{1};
for i=2:length(B)
    B_all=[B_all; B{i}];
end