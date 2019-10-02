function [fig_id]=my_trisurf3(x_full,nodes,elems,map,n_aff,fracture_matrice,node_map_on_double)
rng(0)
x_full=x_full(node_map_on_double);
fig_id=figure(102);

tmp=reshape(x_full,2,length(nodes))';
vel_posunu=sqrt(sum(tmp.^2,2));
nodes=reshape(x_full,2,length(nodes))'+nodes;
p=nodes';
t=elems';
x=p(1,:);
y=p(2,:);
P=[x(t(:));y(t(:))];
T=reshape(1:size(P,2),[3 size(P,2)/3]);
% create random u for testing
if max(map)==0
    tmp=-[map;map;map];
else
    
    tmp=-[map;map;map]/max(map);
end

h=trisurf(elems,x,y,vel_posunu);
h.EdgeColor = 'none';
h.FaceAlpha=0.5;
colormap jet(1000)
axis equal
view(0,90)

hold on

idx_n=sum(n_aff,2)>1;
plot(nodes(idx_n,1),nodes(idx_n,2),'k.','MarkerSize',5)

for i=1:length(fracture_matrice)
    tmp=fracture_matrice{i};
    tmpx=[x(tmp.above_nodes) zeros(length(tmp.above_nodes),1)/0]';
    tmpy=[y(tmp.above_nodes) zeros(length(tmp.above_nodes),1)/0]';
    plot(tmpx(:),tmpy(:),'k-','LineWidth',0.5)
    tmpx=[x(tmp.under_nodes) zeros(length(tmp.above_nodes),1)/0]';
    tmpy=[y(tmp.under_nodes) zeros(length(tmp.above_nodes),1)/0]';
    plot(tmpx(:),tmpy(:),'k-','LineWidth',0.5)
end
title('FETI solution')
hold off
end