function [ fractures ] = duplicate_node_in_other_fractures( fractures,current_frac_id,old_node_id,new_node_id )
%DUPLICATE_NODE_IN_OTHER_FRACTURES Summary of this function goes here
%   Detailed explanation goes here
n=length(fractures);

for i=(current_frac_id+1):n
    fracture=fractures{i};
    mm=length(fracture);
    for j=1:mm
       if(sum(fracture{j}==old_node_id)>0)
           fractures{i}{j}=[fracture{j} new_node_id];
       end
    end
    
end

end

