function [ fracture_matrice ] = update_fracture_matrice( fracture_matrice,node,current_frac_id,node_id,id_new,node_prev,node_next )
%UPDATE_ Summary of this function goes here
%   Detailed explanation goes here
n=length(fracture_matrice);

for i=[1:(current_frac_id-1) (current_frac_id+1):n]
    
    elem_above = fracture_matrice{i}.above_nodes;
    mm=size(elem_above,1);
    for j=1:mm
       if elem_above(j,1)==node_id
           t=(node(elem_above(j,2),:)+node(node_id,:))/2; 
           if mesh.my_elem_classification( node_prev,node(node_id,:),node_next,t )
               fracture_matrice{i}.above_nodes(j,1)=id_new;
           end
       elseif elem_above(j,2)==node_id
           t=(node(elem_above(j,1),:)+node(node_id,:))/2;
           if mesh.my_elem_classification( node_prev,node(node_id,:),node_next,t )
               fracture_matrice{i}.above_nodes(j,2)=id_new;
           end
       end    
    end
    
    elem_under = fracture_matrice{i}.under_nodes;
    for j=1:mm
       if elem_under(j,1)==node_id
           t=(node(elem_under(j,2),:)+node(node_id,:))/2; 
           if mesh.my_elem_classification( node_prev,node(node_id,:),node_next,t )
               fracture_matrice{i}.under_nodes(j,1)=id_new;
           end
       elseif elem_under(j,2)==node_id
           t=(node(elem_under(j,1),:)+node(node_id,:))/2;
           if mesh.my_elem_classification( node_prev,node(node_id,:),node_next,t )
               fracture_matrice{i}.under_nodes(j,2)=id_new;
           end
       end    
    end
    
end

end

