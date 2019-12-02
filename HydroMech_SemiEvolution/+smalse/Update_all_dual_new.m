function [x_elast] = Update_all_dual_new(lambda_ker,idx_no_bounds,idx_bounds,c,B_e,B_i,A_plus,b,R,update_temp_struct)
%UPDATE_B_DUAL Summary of this function goes here
%   Detailed explanation goes here
lambda_ImGt=update_temp_struct.lambda_ImGt;
lambda=lambda_ker+lambda_ImGt;
B=[B_e;-B_i];
x_elast=A_plus*(b-B'*lambda);

% sestaveni alpha jako indexy doplnku x v kernelu A
idx_active_equality=false(length(c),1);
idx_active_equality(idx_no_bounds)=true;
idx_active_equality(idx_bounds)=lambda(idx_bounds)>0; % lambda>0 = kontakt
B_pruh=B(idx_active_equality,:);
c_pruh=c(idx_active_equality);
B_R=B_pruh*R;

alpha=((full(B_R'*B_R)))\(B_R'*(c_pruh-B_pruh*x_elast));
x_elast=x_elast+R*alpha;
end