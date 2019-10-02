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