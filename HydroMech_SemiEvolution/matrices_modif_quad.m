function [ B,G,Au ] = matrices_modif_quad( fracture_matrice,node )
%UPDATE_ Summary of this function goes here
% resulting rectangular matrix with dimensions:
% N ... nodes of 2d domain + duplicities (length of node)
% N_f ... nodes of 1d fracture network + duplicities (length of node_f)

% required data for each fracture side:
% elem_d ... (node_d_1 node_d_2) in domain numbering 
% elem_f ... (node_f_1 node_f_2) in fracture numbering
% ALPHA ... alpha coefficient

N_d_node = size(node,1);
N_f_node = 0;
N_f_elem = 0;
for i=1:length(fracture_matrice)
    N_f_node = N_f_node + length(fracture_matrice{i}.under_material)+1;
    N_f_elem = N_f_elem + (length(fracture_matrice{i}.under_material)*2+1)*2;
end
M_local = [1/6 1/3 0; 0 1/3 1/6];
M2_local = [4 2 -1; 2 16 2; -1 2 4]/30;
B = sparse(N_f_elem, N_d_node);
G = sparse(N_f_node, N_f_elem);
Au = sparse(N_f_elem, N_f_elem);

n=length(fracture_matrice);
counter_f_node=0;
counter_f_elem=0;
for fracture_id=1:n
    for side=1:2
        if side==1
            idx_d_node=fracture_matrice{fracture_id}.under_nodes;
            ALPHA=fracture_matrice{fracture_id}.under_material;
            len=length(ALPHA);
            idx_f_node=[((counter_f_node+1):(counter_f_node+len))' ((counter_f_node+2):(counter_f_node+len+1))'];
            counter_f_node = counter_f_node+len+1;
        else
            idx_d_node=fracture_matrice{fracture_id}.above_nodes;
            ALPHA=fracture_matrice{fracture_id}.above_material;
        end
        idx_f_elem=[((counter_f_elem+1):2:(counter_f_elem+2*len))' ((counter_f_elem+2):2:(counter_f_elem+2*len+1))' ((counter_f_elem+3):2:(counter_f_elem+2*len+2))'];
        counter_f_elem = counter_f_elem+2*len+1;
        LENGTHS = node(idx_d_node(:,2),:) - node(idx_d_node(:,1),:);
        LENGTHS = sqrt(sum(LENGTHS.^2,2));
        for i=1:3
            for j=1:2
                ii = idx_f_elem(:,i);
                jj = idx_d_node(:,j);
                B = B + sparse(ii,jj,M_local(j,i)*LENGTHS,N_f_elem,N_d_node);
    %             B = B + sparse(ii,jj,M_local(j)*LENGTHS.*ALPHA,N_f_elem,N_d_node);
            end
        end
        for i=1:2
            for j=1:3
                ii = idx_f_node(:,i);
                jj = idx_f_elem(:,j);
                G = G + sparse(ii,jj,M_local(i,j)*LENGTHS,N_f_node, N_f_elem);
    %             G = G + sparse(ii,jj,M_local(i)*LENGTHS.*ALPHA,N_f_node, N_f_elem);
            end
        end
        for i=1:3
            for j=1:3
                ii = idx_f_elem(:,i);
                jj = idx_f_elem(:,j);
                Au = Au + sparse(ii,jj,M2_local(i,j)*LENGTHS./ALPHA,N_f_elem,N_f_elem);
            end
        end
    end
end
