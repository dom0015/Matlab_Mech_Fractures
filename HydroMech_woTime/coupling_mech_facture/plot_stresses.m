function [fig_id1,fig_id2,fig_id3] = plot_stresses(x,stress_matrix,sub_elems,Sub_nodes,NODES,fracture_matrice)
%PLOT_STRESSES Summary of this function goes here
%   Detailed explanation goes here
all_stresses=stress_matrix*x;
stres1=all_stresses(1:3:end);
stres2=all_stresses(2:3:end);
stres3=all_stresses(3:3:end);

n=length(sub_elems);

for i=2:n
    sub_elems{i}=sub_elems{i}+max(sub_elems{i-1}(:));
end
elems=cell2mat(sub_elems);
nodes=cell2mat(Sub_nodes);




fig_id1=figure;
p=nodes';
t=elems';
x=p(1,:);
y=p(2,:);
P=[x(t(:));y(t(:))];
T=reshape(1:size(P,2),[3 size(P,2)/3]);
% create random u for testing   
 tmp=[stres1 stres1 stres1]';
h=trisurf(T',P(1,:),P(2,:),tmp(:));
h.EdgeColor = 'none';
h.FaceAlpha=0.5;
colormap jet(1000)
axis equal
view(0,90)
hold on

% idx_n=sum(n_aff,2)>1;
% plot(nodes(idx_n,1),nodes(idx_n,2),'k.','MarkerSize',5)
p=NODES';
x=p(1,:);
y=p(2,:);

for i=1:length(fracture_matrice)
    tmp=fracture_matrice{i};
    tmpx=[x(tmp.above_nodes) zeros(length(tmp.above_nodes),1)/0]';
    tmpy=[y(tmp.above_nodes) zeros(length(tmp.above_nodes),1)/0]';
    plot(tmpx(:),tmpy(:),'k-','LineWidth',0.5)
    tmpx=[x(tmp.under_nodes) zeros(length(tmp.above_nodes),1)/0]';
    tmpy=[y(tmp.under_nodes) zeros(length(tmp.above_nodes),1)/0]';
    plot(tmpx(:),tmpy(:),'k-','LineWidth',0.5)
end
hold off


fig_id2=figure;
p=nodes';
t=elems';
x=p(1,:);
y=p(2,:);
P=[x(t(:));y(t(:))];
T=reshape(1:size(P,2),[3 size(P,2)/3]);
% create random u for testing   
 tmp=[stres2 stres2 stres2]';
h=trisurf(T',P(1,:),P(2,:),tmp(:));
h.EdgeColor = 'none';
h.FaceAlpha=0.5;
colormap jet(1000)
axis equal
view(0,90)
hold on
p=NODES';
x=p(1,:);
y=p(2,:);

for i=1:length(fracture_matrice)
    tmp=fracture_matrice{i};
    tmpx=[x(tmp.above_nodes) zeros(length(tmp.above_nodes),1)/0]';
    tmpy=[y(tmp.above_nodes) zeros(length(tmp.above_nodes),1)/0]';
    plot(tmpx(:),tmpy(:),'k-','LineWidth',0.5)
    tmpx=[x(tmp.under_nodes) zeros(length(tmp.above_nodes),1)/0]';
    tmpy=[y(tmp.under_nodes) zeros(length(tmp.above_nodes),1)/0]';
    plot(tmpx(:),tmpy(:),'k-','LineWidth',0.5)
end
hold off

fig_id3=figure;
p=nodes';
t=elems';
x=p(1,:);
y=p(2,:);
P=[x(t(:));y(t(:))];
T=reshape(1:size(P,2),[3 size(P,2)/3]);
% create random u for testing   
 tmp=[stres3 stres3 stres3]';
h=trisurf(T',P(1,:),P(2,:),tmp(:));
h.EdgeColor = 'none';
h.FaceAlpha=0.5;
colormap jet(1000)
axis equal
view(0,90)
hold on
p=NODES';
x=p(1,:);
y=p(2,:);

for i=1:length(fracture_matrice)
    tmp=fracture_matrice{i};
    tmpx=[x(tmp.above_nodes) zeros(length(tmp.above_nodes),1)/0]';
    tmpy=[y(tmp.above_nodes) zeros(length(tmp.above_nodes),1)/0]';
    plot(tmpx(:),tmpy(:),'k-','LineWidth',0.5)
    tmpx=[x(tmp.under_nodes) zeros(length(tmp.above_nodes),1)/0]';
    tmpy=[y(tmp.under_nodes) zeros(length(tmp.above_nodes),1)/0]';
    plot(tmpx(:),tmpy(:),'k-','LineWidth',0.5)
end
hold off
end

