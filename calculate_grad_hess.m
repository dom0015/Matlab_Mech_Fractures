function [grad_P,hess_P,crossed_L,crossed_R] = calculate_grad_hess(u, POINTS, fracture_matrice,NODE_ABOVE,NODE_UNDER)
%PENALTY_MATRIX Summary of this function goes here
%   Detailed explanation goes here

[OVERLAP_L,OVERLAP_R] = calculate_overlap(u, POINTS, fracture_matrice,NODE_ABOVE,NODE_UNDER);

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
hess_P=sparse(length(u),length(u));
grad_P=sparse(length(u),1);

no_fractures=length(fracture_matrice);

%% normaly, smery otevreni, delky
for f=1:no_fractures
    NODE_above=NODE_ABOVE{f};
    NODE_under=NODE_UNDER{f};
    overlap_L=OVERLAP_L{f};
    overlap_R=OVERLAP_R{f};
    crossed_L=(overlap_L>0);
    crossed_R=(overlap_R>0);

    ua_L=NODE_above(1:end-1,:);
    uu_L=NODE_under(1:end-1,:);
    ua_R=NODE_above(2:end,:);
    uu_R=NODE_under(2:end,:);
    temp=[ua_L uu_L ua_R uu_R];
    indices=[2*temp-1 2*temp];
    
    for i=1:length(ua_L)
        idx=indices(i,:);
        Hd_L_u=M_L*u(idx);
        Hd_R_u=M_R*u(idx);
        temp_L=Hd_L_u*(Hd_L_u');% + overlap_L(i)*M_L;
        temp_R=Hd_R_u*(Hd_R_u');% + overlap_R(i)*M_R;
        temp=2*crossed_L(i)*temp_L + 2*crossed_R(i)*temp_R;
        v=temp(:);
        grad_P=grad_P+sparse(idx,idx*0+1,2*crossed_L(i)*overlap_L(i)*Hd_L_u+...
                                         2*crossed_R(i)*overlap_R(i)*Hd_R_u,length(u),1);
        idx=repmat(idx,length(idx),1);
        jj=idx(:);
        idx=idx';
        ii=idx(:);
        hess_P=hess_P+sparse(ii,jj,v,length(u),length(u));
    end
end
