function [B_all] = calculate_B_equality(u, POINTS, fracture_matrice,NODE_ABOVE,NODE_UNDER,CROSSED_L,CROSSED_R)
%PENALTY_MATRIX Summary of this function goes here
%   Detailed explanation goes here

temp1=u(1:2:end);
temp2=u(2:2:end);
u_reshaped=[temp1 temp2]+POINTS;
u=u_reshaped'; u=u(:);

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
    NODE_above=NODE_ABOVE{f};
    NODE_under=NODE_UNDER{f};
    crossed_L=CROSSED_L{f};
    crossed_R=CROSSED_R{f};

    ua_L=NODE_above(1:end-1,:);
    uu_L=NODE_under(1:end-1,:);
    ua_R=NODE_above(2:end,:);
    uu_R=NODE_under(2:end,:);
    temp=[ua_L uu_L ua_R uu_R];
    indices=[2*temp-1 2*temp];

    N=sum(crossed_L)+sum(crossed_R);
    mat_B=sparse(N,length(u));
    i_B=1;

    for i=1:length(ua_L)
        idx=indices(i,:);
        Hd_L_u=M_L*u(idx);
        Hd_R_u=M_R*u(idx);
        if crossed_L(i) % matrix B for equality-constrained optimization
            mat_B(i_B,idx)=Hd_L_u; i_B=i_B+1;
        end
        if crossed_R(i)
            mat_B(i_B,idx)=Hd_R_u; i_B=i_B+1;
        end
    end
    B{f}=mat_B;
end
B_all=B{1};
for i=2:length(B)
    B_all=[B_all; B{i}];
end