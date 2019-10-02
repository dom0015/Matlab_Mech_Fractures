function [] = plot_divergence(u,node,elem)
%PLOT_GRAD Summary of this function goes here
%   Detailed explanation goes here
p=node';
t=elem';
x=p(1,:);
y=p(2,:);
P=[x(t(:));y(t(:))];
T=reshape(1:size(P,2),[3 size(P,2)/3]);
% create random u for testing
tmp=[u';u';u'];
h=trisurf(T',P(1,:),P(2,:),tmp(:));
h.EdgeColor = 'none';
colormap jet(1000)
axis equal
view(0,90)
colorbar