function [fig_id]=my_trisurf2(x_full,nodes,elems,map,n_aff,fracture_matrice,node_map_on_double)
rng(0)
x_full=x_full(node_map_on_double);
fig_id=figure;
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

h=trisurf(T',P(1,:),P(2,:),tmp(:));
h.EdgeColor = 'none';
h.FaceAlpha=0.5;
colormap HSV(1000)
colors=max(map)+1;
t_s=[0.7 1];
S=t_s(mod(0:colors-1,2)+1);
t_v=[0.85 0.85 1 1];
V=t_v(mod(0:colors-1,4)+1);
colors_corr=ceil(colors/1);
t_h=linspace(0,1-1/colors_corr,colors_corr);
t_h=t_h(floor((0:colors-1)/1)+1);
[~,idx]=sort(rand(colors,1));
H=t_h(idx);
V=V(idx);
S=S(idx);
col=hsv2rgb([H' S' V']);
colormap(col)
axis equal
view(0,90)

hold on

idx_n=sum(n_aff,2)>1;
plot(nodes(idx_n,1),nodes(idx_n,2),'k.','MarkerSize',3)

for i=1:length(fracture_matrice)
    tmp=fracture_matrice{i};
    tmpx=[x(tmp.above_nodes) zeros(length(tmp.above_nodes),1)/0]';
    tmpy=[y(tmp.above_nodes) zeros(length(tmp.above_nodes),1)/0]';
    plot(tmpx(:),tmpy(:),'k-','LineWidth',0.5)
    tmpx=[x(tmp.under_nodes) zeros(length(tmp.above_nodes),1)/0]';
    tmpy=[y(tmp.under_nodes) zeros(length(tmp.above_nodes),1)/0]';
    plot(tmpx(:),tmpy(:),'k-','LineWidth',0.5)
end
end