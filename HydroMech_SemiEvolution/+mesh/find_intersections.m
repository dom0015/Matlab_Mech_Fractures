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