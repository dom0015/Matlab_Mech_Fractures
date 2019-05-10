function [ fractures ] = remap_fracture_nodeid( fractures,node_map )
%REMAP_FRACTURE_NODEID Summary of this function goes here
%   Detailed explanation goes here
n=length(fractures);

for i=1:n
    mm=length(fractures{i});
    for j=1:mm
        fractures{i}{j}=node_map(fractures{i}{j});
    end
    
end


end

