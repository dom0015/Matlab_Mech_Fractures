function [grad_P,hess_P,NODE_above,NODE_under,crossed_L,crossed_R,mat_B] = penalty_matrix_exact_approximated(u, POINTS, fracture_matrice,intersections)
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
NODE_above = [NODE_start; NODE_above; NODE_end];
NODE_under = [NODE_start; NODE_under; NODE_end];
COORD_above=u_reshaped(NODE_above,:);
COORD_under=u_reshaped(NODE_under,:);
COORD_middles=(COORD_under+COORD_above)/2;
%COORD_middles=POINTS([NODE_start; NODE_above; NODE_end],:);

temp=(COORD_middles(2:end,:)-COORD_middles(1:end-1,:));
NORMALS_scaled=[temp(:,2) -temp(:,1)];
LENGTHS=sqrt(sum(temp.^2,2));
NORMALS=NORMALS_scaled./repmat(LENGTHS,1,2);
DISPLACEMENTS=COORD_above-COORD_under;

overlap_L=dot(NORMALS_scaled,DISPLACEMENTS(1:end-1,:),2);
overlap_R=dot(NORMALS_scaled,DISPLACEMENTS(2:end,:),2);
crossed_L=(overlap_L>0);
crossed_R=(overlap_R>0);

%figure; plot(overlap_L); hold on; plot(overlap_L.*crossed_L); plot(overlap_R); plot(overlap_R.*crossed_R); 

Z=zeros(4);
M_R_variable=[0 0 -0.5 0.5; 0 0 -0.5 0.5; 0   0    0.5 -0.5;  0    0   0.5 -0.5];
M_R_constant=[0 0  0   0  ; 0 0  0   0  ; 0.5 0.5 -0.5 -0.5; -0.5 -0.5 0.5  0.5];
M_L_variable=[-0.5 0.5  0    0  ; -0.5  0.5 0   0  ; 0.5 -0.5 0 0; 0.5 -0.5 0 0];
M_L_constant=[ 0.5 0.5 -0.5 -0.5; -0.5 -0.5 0.5 0.5; 0    0   0 0; 0    0   0 0];
M_R=M_R_constant;%+M_R_variable;
M_L=M_L_constant;%+M_L_variable;
M_R=[Z -M_R; M_R Z];
M_L=[Z -M_L; M_L Z];
hess_P=sparse(length(u),length(u));
grad_P=sparse(length(u),1);

ua_L=NODE_above(1:end-1,:);
uu_L=NODE_under(1:end-1,:);
ua_R=NODE_above(2:end,:);
uu_R=NODE_under(2:end,:);
temp=[ua_L uu_L ua_R uu_R];
%temp=[ua_R uu_R ua_L uu_L];
%temp=[uu_L ua_L uu_R ua_R];
indices=[2*temp-1 2*temp];
indicesH=[uu_L ua_L uu_R ua_R];

N=sum(crossed_L)+sum(crossed_R);
mat_B=sparse(N,length(u));
i_B=1;

for i=1:length(LENGTHS)
    l=LENGTHS(i);
    idx=indices(i,:);
    Hd_L_u=M_L*u(idx);
    Hd_R_u=M_R*u(idx);
    if crossed_L(i) % matrix B for equality-constrained optimization
        mat_B(i_B,idx)=Hd_L_u; i_B=i_B+1;
    end
    if crossed_R(i)
        mat_B(i_B,idx)=Hd_R_u; i_B=i_B+1;
    end
    temp_L=Hd_L_u*(Hd_L_u');% + overlap_L(i)*M_L;
    temp_R=Hd_R_u*(Hd_R_u');% + overlap_R(i)*M_R;
%     temp_L=Hd_L_u*(Hd_L_u') + M_L;
%     temp_R=Hd_R_u*(Hd_R_u') + M_R;
%     temp_L=Hd_L_u*(Hd_L_u');
%     temp_R=Hd_R_u*(Hd_R_u');
    temp=2*crossed_L(i)*temp_L + 2*crossed_R(i)*temp_R;
    v=temp(:);
    grad_P=grad_P+sparse(idx,idx*0+1,2*crossed_L(i)*overlap_L(i)*Hd_L_u+...
                                     2*crossed_R(i)*overlap_R(i)*Hd_R_u,length(u),1);
%     grad_P=grad_P+sparse(idx,idx*0+1,2*crossed_L(i)*(bla)*Hd_L_u+...
%                                      2*crossed_R(i)*(bla)*Hd_R_u,length(u),1);
    idx=repmat(idx,length(idx),1);
    %idx=repmat(iH(:)',length(idx),1);
    %v=2*l*crossed_L(i)*temp_L_H(:);
    jj=idx(:);
    idx=idx';
    ii=idx(:);
    hess_P=hess_P+sparse(ii,jj,v,length(u),length(u));
end

end

