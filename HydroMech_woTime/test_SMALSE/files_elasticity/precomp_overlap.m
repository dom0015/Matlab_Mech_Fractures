function [NODE_ABOVE,NODE_UNDER] = precomp_overlap(fracture_matrice,intersections)
%CALCULATE_OVERLAP Summary of this function goes here
%   Detailed explanation goes here

no_fractures=length(fracture_matrice);

NODE_ABOVE=cell(no_fractures,1);
NODE_UNDER=cell(no_fractures,1);
%% normaly, smery otevreni, delky
for f=1:no_fractures
    ELEM_above=fracture_matrice{f}.above_nodes;
    ELEM_under=fracture_matrice{f}.under_nodes;
    NODE_start=ELEM_above(1,1);
    NODE_end=ELEM_above(end,2);
    NODE_above=ELEM_above(2:end,1);
    NODE_under=ELEM_under(2:end,1);
    if ~isempty(intersections)
        NODE_above_2=ELEM_above(1:end-1,2);
        temp_a=(NODE_above~=NODE_above_2); % temp is nonzero in case of intersecting fractures
        NODE_under_2=ELEM_under(1:end-1,2);
%         temp_u=(NODE_under~=NODE_under_2); % temp is nonzero in case of intersecting fractures
        shift=0;
        for i=find(temp_a)'
            NODE_above=[NODE_above(1:(i+shift-1)); NODE_above_2(i); NODE_above((i+shift):end)];
            NODE_under=[NODE_under(1:(i+shift-1)); NODE_under_2(i); NODE_under((i+shift):end)];
            shift=shift+1;
        end
    end
    NODE_ABOVE{f} = [NODE_start; NODE_above; NODE_end];
    NODE_UNDER{f} = [NODE_start; NODE_under; NODE_end];
end