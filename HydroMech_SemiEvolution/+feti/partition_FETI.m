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
    [map,edgecut] = metis.metismex('PartGraphRecursive',W,partitions,opts);
else
    edgecut=0;
    map=zeros(length(elems),1);
end
fprintf('edgecut=%d\n',edgecut);

[sub_nodes,sub_elem,sub_material_constants,sub_volume_force,sub_boundary,...
    sub_boundary_val,sub_sizes,glueB,node_map_on_double,n_aff]=...
    feti.subdomains_geometry(nodes,elems,material_constants,volume_force,...
    Neumann_boundaries,N_bound_value,map);

divide_neumann_boundary=@(N_bound,N_bound_value)feti.subdomains_Neumann(N_bound,N_bound_value,map);
divide_neumann_boundary_and_force=@(N_bound,N_bound_value,volume_force)feti.subdomains_Neumann_force(N_bound,N_bound_value,volume_force,map);

%% plotting ---------------------------------------------------------------
plot_func=@(x,y)viz.my_trisurf2(x,nodes,elems,map,n_aff,y,node_map_on_double);
plot_func2=@(x,y)viz.my_trisurf3(x,nodes,elems,map,n_aff,y,node_map_on_double);
end
