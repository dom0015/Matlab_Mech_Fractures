function [ first ] = my_elem_classification( ps,pm,pe,t )
%MY_ELEM_CLASSIFICATION Summary of this function goes here
%   Detailed explanation goes here

n=size(t,1);

ps=ps-pm;
pe=pe-pm;
points=t-repmat(pm,n,1);


a_phi=mesh.my_atan(ps(1),ps(2));
b_phi=mesh.my_atan(pe(1),pe(2));

p_phi=mesh.my_atan(points(:,1),points(:,2));

b_corr=mesh.my_angle_diff(b_phi,a_phi);
p_corr=mesh.my_angle_diff(p_phi,a_phi);

first=p_corr<b_corr;

end