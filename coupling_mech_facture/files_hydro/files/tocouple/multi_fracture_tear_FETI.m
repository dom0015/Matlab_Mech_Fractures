function [ node,elem,bdflag,fractures,fracture_matrice,Db] = multi_fracture_tear( node,elem,fractures,bdflag,fracture_matrice,Db)
%ONE_FRACTURE_TEAR Summary of this function goes here
%   Detailed explanation goes here
m=length(fractures);
elem_t(:,1)=node(elem(:,1),1)+node(elem(:,2),1)+node(elem(:,3),1);
elem_t(:,2)=node(elem(:,1),2)+node(elem(:,2),2)+node(elem(:,3),2);
elem_t=elem_t/3;

for fracture_id=1:m
    fracture=fractures{fracture_id};
    n=length(fracture);
    nodes_count=length(node);
    nodes_map=1:nodes_count;
    
    for i=2:n-1
        nn=length(fracture{i});
        
        for nd=1:nn
            node_id=fracture{i}(nd);
            cor_elems=sum(elem==node_id,2)>0;
            [ under_fracture ] = my_elem_classification( node(fracture{i-1}(1),:),node(node_id,:),node(fracture{i+1}(1),:),elem_t(cor_elems,:));
            if (sum(under_fracture)==0)||(sum(under_fracture)==length(under_fracture))
                
            else
            [ fractures ] = duplicate_node_in_other_fractures( fractures,fracture_id,node_id,nodes_count+1 );
            %soucasna trhlina under precislovat na node_id na nodes_count+1
            fracture_matrice{fracture_id}.under_nodes(fracture_matrice{fracture_id}.under_nodes==node_id)=nodes_count+1; 
            %vsechny ostatni trhliny projit a tam kde se cokoli rovna
            %node_id tak pomoci duplicate_node_in_other_fractures
            fracture_matrice = update_fracture_matrice(fracture_matrice,node,fracture_id,node_id,nodes_count+1,node(fracture{i-1}(1),:),node(fracture{i+1}(1),:));

            node=[node;node(node_id,:)];
            cor_elems(cor_elems)=under_fracture;
            elem(repmat(cor_elems,1,3)&(elem==node_id))=nodes_count+1;
            node_id=nodes_map(node_id);
            nodes_map(nodes_map>node_id)=nodes_map(nodes_map>node_id)+1;
            nodes_map=[nodes_map node_id+1];
            nodes_count=nodes_count+1;
            end
        end
    end
    node(nodes_map,:)=node;
    elem=nodes_map(elem);
    

    [ fractures ] = remap_fracture_nodeid( fractures,nodes_map );
    %precislovat nody na trhlinach dle mapy
    for fid=1:m
        fracture_matrice{fid}.under_nodes=nodes_map(fracture_matrice{fid}.under_nodes);
        fracture_matrice{fid}.above_nodes=nodes_map(fracture_matrice{fid}.above_nodes);
    end
    
    if nargin==6
        Db=[Db false(4,length(nodes_map)-length(Db))];
        Db(:,nodes_map)=Db;
    end
end

end

