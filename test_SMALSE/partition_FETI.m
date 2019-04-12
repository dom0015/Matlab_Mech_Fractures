function [map,sub_nodes,sub_elem,sub_material_constants,sub_volume_force,sub_boundary,sub_boundary_val,sub_sizes,glueB,node_map_on_double]=partition_FETI(nodes,elems,material_constants,volume_force,Neumann_boundaries,N_bound_value,partitions)

n=length(nodes);
m=length(elems);

B=sparse([1:m 1:m 1:m],[elems(:,1) elems(:,2) elems(:,3)],ones(3*m,1),m,n);
W=(B*B')>=2;
opts.seed=13;
opts.niter=10;
opts.ncuts=2;
[map,edgecut] = metismex('PartGraphRecursive',W,partitions,opts);
fprintf('edgecut=%d\n',edgecut);
[n_aff]=node_affinity(nodes,elems,map);

[sub_nodes,sub_elem,sub_material_constants,sub_volume_force,sub_boundary,sub_boundary_val,sub_sizes,glueB,node_map_on_double]=subdomains_geometry(nodes,elems,material_constants,volume_force,Neumann_boundaries,N_bound_value,map,n_aff);

%% plotting ---------------------------------------------------------------
my_trisurf(nodes,elems,map);
hold on
idx_n=sum(n_aff,2)>1;
plot(nodes(idx_n,1),nodes(idx_n,2),'k.','MarkerSize',10)
plot(nodes(idx_n,1),nodes(idx_n,2),'w.','MarkerSize',1)
end

%%
function [n_aff]=node_affinity(nodes,elems,map)
n=length(nodes);
m=max(map)+1;
n_aff=false(n,m);
for i=1:m
    locnodes=elems(map==(i-1),:);
    n_aff(locnodes,i)=true;
end
end

%%
function [sub_nodes,sub_elem,sub_material_constants,sub_volume_force,sub_boundary,sub_boundary_val,sub_sizes,glueB_double,node_map_on_double]=subdomains_geometry(nodes,elems,material_constants,volume_force,Neumann_boundaries,N_bound_value,map,n_aff)
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
n_new_nodes=sum(n_aff(:));
node_map_no=zeros(n_new_nodes,1);
node_map_on=zeros(length(nodes),1);
for i=1:n
    sub_sizes(i,1)=sum(n_aff(:,i));
    if i==1
        sub_sizes(i,2)=1;
        sub_sizes(i,3)=sub_sizes(i,1);
    else
        sub_sizes(i,2)=sub_sizes(i-1,3)+1;
        sub_sizes(i,3)=sub_sizes(i-1,3)+sub_sizes(i,1);
    end
    sub_nodes{i}=nodes(n_aff(:,i),:);
    tmp=elems(map==(i-1),:);
    tmp_idx_holder(n_aff(:,i))=1:sub_sizes(i,1);
    sub_elem{i}=tmp_idx_holder(tmp);
    sub_boundary{i}=Neumann_boundaries(map==(i-1),:);
    tmp_cell{1}=N_bound_value{1}(map==(i-1),:);
    tmp_cell{2}=N_bound_value{2}(map==(i-1),:);
    sub_boundary_val{i}=tmp_cell;
    sub_material_constants{i}=material_constants(map==(i-1),:);
    sub_volume_force{i}=volume_force(map==(i-1),:);
    tmp_loc_nodes_idx{i}=sub_sizes(i,2):sub_sizes(i,3);
    node_map_no(sub_sizes(i,2):sub_sizes(i,3))=n_idx(n_aff(:,i));
    node_map_on(n_idx(n_aff(:,i)))=sub_sizes(i,2):sub_sizes(i,3);
end
nn=sub_sizes(end,3);
nn_aff=sum(n_aff,2);
mm=max(nn_aff);
tmp_sum=0;

for i=2:mm
    tmp_sum=tmp_sum+i*(i-1)/2*sum(nn_aff==i);
end


% temp_glue_j=zeros(tmp_sum,2);
% temp_nidx=zeros(tmp_sum,1);
% it_idx=0;
% tmp=n_aff(:,1);
% for i=1:n
%     tmp(:)=false;
%     tmp_i_idx=find(n_aff(:,i));
%     for j=1:(i-1)
%         tmp(tmp_i_idx)=n_aff(tmp_i_idx,j);
%         if sum(tmp(tmp_i_idx))>0
%             i_idxs=tmp_loc_nodes_idx{i}(tmp(n_aff(:,i)));
%             j_idxs=tmp_loc_nodes_idx{j}(tmp(n_aff(:,j)));
%             mm=length(i_idxs);
%             temp_nidx((it_idx+1):(it_idx+mm))=n_idx(tmp);
%             temp_glue_j((it_idx+1):(it_idx+mm),:)=[i_idxs' j_idxs'];
%             it_idx=it_idx+mm;
%         end
%     end
% end

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

%glueB=sparse([(1:n)';(1:n)'], temp_glue_j(:),[ones(n,1);-ones(n,1)],n,nn);

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
function [fig_id]=my_trisurf(nodes,elems,map)

fig_id=figure;
p=nodes';
t=elems';
x=p(1,:);
y=p(2,:);
P=[x(t(:));y(t(:))];
T=reshape(1:size(P,2),[3 size(P,2)/3]);
% create random u for testing
tmp=-[map;map;map]/max(map);
h=trisurf(T',P(1,:),P(2,:),tmp(:));
h.EdgeColor = 'none';
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
end