function [ fractures_cells,fracture_matrice,intersections,lengths] = fracture2cells_geometry( fractures )
%FRACTURE2CELLS Summary of this function goes here
%   Detailed explanation goes here
[intersections]=find_intersections(fractures);
n=length(fractures);
lengths=zeros(1,n);
fractures_cells=cell(n,1);
fracture_matrice=cell(n,1);
for i=1:n
    tmp=fractures{i};
    m=length(tmp);
    lengths(i)=m;
    tmp_fm.above_nodes=[tmp(1:end-1)' tmp(2:end)'];
    tmp_fm.under_nodes=[tmp(1:end-1)' tmp(2:end)'];
    fracture_matrice{i}=tmp_fm;
    tmp2=cell(m,1);
    for j=1:m
       tmp2{j}=tmp(j); 
    end
    fractures_cells{i}=tmp2;
end

end

function [intersections]=find_intersections(fractures)
n=length(fractures);
all_nodes=cell2mat(fractures');
[~, ind] = unique(all_nodes);
% duplicate indices
duplicate_ind = setdiff(1:length(all_nodes), ind);
% duplicate values
duplicate_value = unique(all_nodes(duplicate_ind));
intersects=length(duplicate_value);
intersections=zeros(intersects,n);

for i=1:intersects
    for j=1:n
        tmp=find(fractures{j}==duplicate_value(i));
        if isempty(tmp)
           tmp=-1; 
        end
        intersections(i,j)=tmp;
    end
end

end