function [ fractures, fractures_positions, no_fractures ] = create_fractures( frac_start_end, node, h_elem )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
no_fractures=size(frac_start_end,1);
fractures = cell(no_fractures,1);
fractures_positions = cell(no_fractures,1);
for i=1:no_fractures
    fractures{i} = one_fracture_nodes( h_elem, frac_start_end{i,1}, frac_start_end{i,2} );
    fractures_positions{i} = node(fractures{i},:);
end

end

