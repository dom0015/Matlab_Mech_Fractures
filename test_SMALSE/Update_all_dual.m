function [F,ff,G,c_ker,update_temp_struct,geom_res,x_elast] = Update_all_dual(lambda_ker,idx_no_bounds,idx_bounds,c,B_e,A_plus,b,R,B_iupdate,update_temp_struct)
%UPDATE_B_DUAL Summary of this function goes here
%   Detailed explanation goes here
lambda_ImGt=update_temp_struct.lambda_ImGt;
B_i=update_temp_struct.B_i;
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

[~,~,tv]=find(update_temp_struct.B_i-B_i);
geom_res=norm(tv);
update_temp_struct.B_i=B_i;
B=[B_e;-B_i];

F11=update_temp_struct.Be_Ap_Be;
F12=update_temp_struct.Be_Ap*B_i';
F22=B_i*A_plus*B_i';

%F=B*A_plus*B';
F=[F11 -F12;-F12' F22];


d=B*(update_temp_struct.Ap_b);
G=R'*B';
e=R'*b;

lambda_ImGt=G'*((G*G')\(e));

update_temp_struct.lambda_ImGt=lambda_ImGt;

c_ker=-lambda_ImGt;
c_ker(idx_no_bounds)=-Inf;
ff=d-F*lambda_ImGt;
% global tmprint
% tmprint(x_elast);
end