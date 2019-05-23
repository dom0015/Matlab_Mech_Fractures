function [map,sub_nodes,sub_elem,sub_material_constants,sub_volume_force,...
    sub_boundary,sub_boundary_val,sub_sizes,glueB,node_map_on_double,...
    plot_func,plot_func2,divide_neumann_boundary,divide_neumann_boundary_and_force]=partition_FETI...
    (nodes,elems,material_constants,volume_force,Neumann_boundaries,N_bound_value,partitions)

n=length(nodes);
m=length(elems);

B=sparse([1:m 1:m 1:m],[elems(:,1) elems(:,2) elems(:,3)],ones(3*m,1),m,n);
W=(B*B')>=2;
opts.seed=13;
opts.niter=10;
opts.ncuts=2;
if partitions>1
    [map,edgecut] = metismex('PartGraphRecursive',W,partitions,opts);
else
    edgecut=0;
    map=zeros(length(elems),1);
end
fprintf('edgecut=%d\n',edgecut);

[sub_nodes,sub_elem,sub_material_constants,sub_volume_force,sub_boundary,...
    sub_boundary_val,sub_sizes,glueB,node_map_on_double,n_aff]=...
    subdomains_geometry(nodes,elems,material_constants,volume_force,...
    Neumann_boundaries,N_bound_value,map);

divide_neumann_boundary=@(N_bound,N_bound_value)subdomains_Neumann(N_bound,N_bound_value,map);
divide_neumann_boundary_and_force=@(N_bound,N_bound_value,volume_force)subdomains_Neumann_force(N_bound,N_bound_value,volume_force,map);

%% plotting ---------------------------------------------------------------
% my_trisurf(nodes,elems,map);
% hold on
% idx_n=sum(n_aff,2)>1;
% plot(nodes(idx_n,1),nodes(idx_n,2),'k.','MarkerSize',10)
% plot(nodes(idx_n,1),nodes(idx_n,2),'w.','MarkerSize',1)

plot_func=@(x,y)my_trisurf2(x,nodes,elems,map,n_aff,y,node_map_on_double);
plot_func2=@(x,y)my_trisurf3(x,nodes,elems,map,n_aff,y,node_map_on_double);
end

%%
function [sub_nodes,sub_elem,sub_material_constants,sub_volume_force,...
    sub_boundary,sub_boundary_val,sub_sizes,glueB_double,node_map_on_double,n_aff]=...
    subdomains_geometry(nodes,elems,material_constants,volume_force,...
    Neumann_boundaries,N_bound_value,map)
n=max(map)+1;
sub_nodes=cell(n,1);
tmp_loc_nodes_idx=cell(n,1);
sub_boundary=cell(n,1);
sub_boundary_val=cell(n,1);
sub_material_constants=cell(n,1);
sub_volume_force=cell(n,1);
sub_elem=cell(n,1);
sub_sizes=zeros(n,3);
n_idx=1:length(nodes);
tmp_idx_holder=zeros(length(nodes),1);
node_map_on=zeros(length(nodes),1);
n_aff=false(length(nodes),n);


for i=1:n
    tmp_loc_elem=map==(i-1);
    locnodes=elems(tmp_loc_elem,:);
    n_aff(locnodes,i)=true;
    
    sub_sizes(i,1)=sum(n_aff(:,i));
    if i==1
        sub_sizes(i,2)=1;
        sub_sizes(i,3)=sub_sizes(i,1);
    else
        sub_sizes(i,2)=sub_sizes(i-1,3)+1;
        sub_sizes(i,3)=sub_sizes(i-1,3)+sub_sizes(i,1);
    end
    sub_nodes{i}=nodes(n_aff(:,i),:);
    tmp_idx_holder(n_aff(:,i))=1:sub_sizes(i,1);
    sub_elem{i}=tmp_idx_holder(locnodes);
    sub_boundary{i}=Neumann_boundaries(tmp_loc_elem,:);
    tmp_cell{1}=N_bound_value{1}(tmp_loc_elem,:);
    tmp_cell{2}=N_bound_value{2}(tmp_loc_elem,:);
    sub_boundary_val{i}=tmp_cell;
    sub_material_constants{i}=material_constants(tmp_loc_elem,:);
    sub_volume_force{i}=volume_force(tmp_loc_elem,:);
    tmp_loc_nodes_idx{i}=sub_sizes(i,2):sub_sizes(i,3);
    node_map_on(n_idx(n_aff(:,i)))=sub_sizes(i,2):sub_sizes(i,3);
end

n_new_nodes=sum(n_aff(:));
node_map_no=zeros(n_new_nodes,1);
for i=1:n
    node_map_no(sub_sizes(i,2):sub_sizes(i,3))=n_idx(n_aff(:,i));
end

nn=sub_sizes(end,3);
nn_aff=sum(n_aff,2);
mm=max(nn_aff);
tmp_sum=0;

for i=2:mm
    tmp_sum=tmp_sum+i*(i-1)/2*sum(nn_aff==i);
end

new_old_incidency=sparse((1:n_new_nodes)',node_map_no,ones(n_new_nodes,1),n_new_nodes,length(nodes));
[idx_i_no,idx_j_no,~]=find(triu(new_old_incidency*new_old_incidency',1));
temp_glue_j=[idx_j_no idx_i_no];
temp_nidx= node_map_no(idx_j_no);

idx_corners=find(nn_aff>2);
idx_stay=true(tmp_sum,1);
for i=idx_corners'
    entries=(temp_nidx==i);
    tmp=temp_glue_j(entries,:);
    pivot=tmp(1,1);
    idx_stay(entries)=((tmp(:,1)==pivot) | (tmp(:,2)==pivot));
end

temp_nidx=temp_nidx(idx_stay);
temp_glue_j=temp_glue_j(idx_stay,:);
n=length(temp_nidx);

idx_1=[temp_glue_j(:,1)'*2-1;temp_glue_j(:,1)'*2];
idx_2=[temp_glue_j(:,2)'*2-1;temp_glue_j(:,2)'*2];
glueB_double=sparse([(1:2*n)';(1:2*n)'], [idx_1(:);idx_2(:)],[ones(2*n,1);-ones(2*n,1)],2*n,2*nn);
node_map_on_double=[node_map_on'*2-1;node_map_on'*2];
node_map_on_double=node_map_on_double(:);
sub_sizes(:,1)=2*sub_sizes(:,1);
sub_sizes(:,2)=2*sub_sizes(:,2)-1;
sub_sizes(:,3)=2*sub_sizes(:,3);
end

%%
function [sub_boundary,sub_boundary_val]=subdomains_Neumann(Neumann_boundaries,N_bound_value,map)
n=max(map)+1;
sub_boundary=cell(n,1);
sub_boundary_val=cell(n,1);
for i=1:n
    tmp_loc_elem=map==(i-1);
    sub_boundary{i}=Neumann_boundaries(tmp_loc_elem,:);
    tmp_cell{1}=N_bound_value{1}(tmp_loc_elem,:);
    tmp_cell{2}=N_bound_value{2}(tmp_loc_elem,:);
    sub_boundary_val{i}=tmp_cell;
end
end
function [sub_boundary,sub_boundary_val,sub_volume_force]=subdomains_Neumann_force(Neumann_boundaries,N_bound_value,volume_force,map)
n=max(map)+1;
sub_boundary=cell(n,1);
sub_volume_force=cell(n,1);
sub_boundary_val=cell(n,1);
for i=1:n
    tmp_loc_elem=map==(i-1);
    sub_boundary{i}=Neumann_boundaries(tmp_loc_elem,:);
    tmp_cell{1}=N_bound_value{1}(tmp_loc_elem,:);
    tmp_cell{2}=N_bound_value{2}(tmp_loc_elem,:);
    sub_boundary_val{i}=tmp_cell;
    sub_volume_force{i}=volume_force(tmp_loc_elem,:);
end
end


%%
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

%%
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
hold off
end