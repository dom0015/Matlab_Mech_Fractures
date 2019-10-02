function [ node,elem,bdflag,fractures_cell,fracture_matrice,Db,Neumann_normalx,Neumann_normaly,fractures] = smooth_fracture( node,elem,fractures_cell,bdflag,fracture_matrice,Db,Neumann_normalx,Neumann_normaly,fractures)
no_elem = size(elem,1);
idx_new_node = size(node,1) + 1;
idx_new_elem = no_elem + 1;
ii = [1:no_elem; 1:no_elem; 1:no_elem]; ii = ii(:);
jj = elem'; jj = jj(:);
map = sparse(ii,jj,1);
no_fractures = length(fracture_matrice);
above_middles_cell = cell(no_fractures);
for i=1:no_fractures
    under_nodes = fracture_matrice{i}.under_nodes;
    under_middles = [];
    above_nodes = fracture_matrice{i}.above_nodes;
    above_middles = [];
    fracture_length = size(under_nodes,1);
    above_middles_cell{i}=cell(1,1);
    for j=1:fracture_length
        under_side = under_nodes(j,:);
        under_triangle_idx = find((map(:,under_side(1)) + map(:,under_side(2))) == 2);
        under_triangle = elem(under_triangle_idx,:);
        elem(under_triangle_idx,under_triangle==under_side(1)) = idx_new_node;
        under_triangle(under_triangle==under_side(2)) = idx_new_node;
        elem(idx_new_elem,:) = under_triangle;
        node(idx_new_node,:) = mean(node(under_side,:));
        under_middles = [under_middles; idx_new_node];
        idx_new_node = idx_new_node + 1;
        idx_new_elem = idx_new_elem + 1;
        
        above_side = above_nodes(j,:);
        above_triangle_idx = find((map(:,above_side(1)) + map(:,above_side(2))) == 2);
        above_triagle = elem(above_triangle_idx,:);
        elem(above_triangle_idx,above_triagle==above_side(1)) = idx_new_node;
        above_triagle(above_triagle==above_side(2)) = idx_new_node;
        elem(idx_new_elem,:) = above_triagle;
        node(idx_new_node,:) = mean(node(above_side,:));
        above_middles = [above_middles; idx_new_node];
        above_middles_cell{i}{j} = idx_new_node;
        idx_new_node = idx_new_node + 1;
        idx_new_elem = idx_new_elem + 1;
    end
    under_nodes_begin = [under_nodes(:,1)'; under_middles']; under_nodes_begin = under_nodes_begin(:);
    under_nodes_end   = [under_middles'; under_nodes(:,2)']; under_nodes_end   = under_nodes_end(:);
    fracture_matrice{i}.under_nodes = [under_nodes_begin under_nodes_end];
    
    above_nodes_begin = [above_nodes(:,1)'; above_middles']; above_nodes_begin = above_nodes_begin(:);
    above_nodes_end   = [above_middles'; above_nodes(:,2)']; above_nodes_end   = above_nodes_end(:);
    fracture_matrice{i}.above_nodes = [above_nodes_begin above_nodes_end];
    
    temp = fractures{i};
    temp2 = [temp(1:end-1); above_middles']; temp2=temp2(:);
    temp = [temp2; temp(end)];
    fractures{i} = temp';
    
    temp = fractures_cell{i}(1:end-1);
    temp2 = [temp'; above_middles_cell{1}]; temp2=temp2(:);
    temp2 = [temp2; fractures_cell{i}(end)];
    fractures_cell{i} = temp2;
    
end

Db=[Db false(4,length(node)-length(Db))];
bdflag=[bdflag; zeros(length(elem)-length(bdflag),3)];
Neumann_normalx=[Neumann_normalx; zeros(length(elem)-length(Neumann_normalx),3)];
Neumann_normaly=[Neumann_normaly; zeros(length(elem)-length(Neumann_normaly),3)];
