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