function [F,ff,G,c_ker] = Update_B_dual(lambda_ker,idx_no_bounds,idx_bounds,c,B_e,A_plus,b,R,B_iupdate)
%UPDATE_B_DUAL Summary of this function goes here
%   Detailed explanation goes here
global lambda_ImGt B_i
lambda=lambda_ker+lambda_ImGt;
B=[B_e;-B_i];
x_elast=A_plus*(b-B'*lambda);

% sestaveni alpha jako indexy doplnku x v kernelu A
idx_active_equality=false(length(c),1);
idx_active_equality(idx_no_bounds)=true;
idx_active_equality(idx_bounds)=lambda(idx_bounds)>0;
B_pruh=B(idx_active_equality,:);
c_pruh=c(idx_active_equality);
B_R=B_pruh*R;
alpha=(B_R'*B_R)\(B_R'*(c_pruh-B_pruh*x_elast));

x_elast=x_elast+R*alpha;

B_i=B_iupdate(x_elast);

B=[B_e;-B_i];
F=B*A_plus*B';
d=B*A_plus*b;
G=R'*B';
e=R'*b;

lambda_ImGt=G'*((G*G')\(e));
c_ker=-lambda_ImGt;
c_ker(idx_no_bounds)=-Inf;
ff=d-F*lambda_ImGt;
end