function [ A,M ] = interaction_matrix( fracture_matrice,node )
%UPDATE_ Summary of this function goes here
% resulting rectangular matrix with dimensions:
% N ... nodes of 2d domain + duplicities (length of node)
% N_f ... nodes of 1d fracture network + duplicities (length of node_f)

% required data for each fracture side:
% elem_d ... (node_d_1 node_d_2) in domain numbering 
% elem_f ... (node_f_1 node_f_2) in fracture numbering
% ALPHA ... alpha coefficient

N = size(node,1);
N_f = 0;
for i=1:length(fracture_matrice)
    N_f = N_f + length(fracture_matrice{i}.under_material)+1;
end
M_local = [1/3 1/6; 1/6 1/3];
A = sparse(N_f, N);
M = sparse(N, N);

n=length(fracture_matrice);
counter=0;
for fracture_id=1:n
    for side=1:2
        if side==1
            elem_d=fracture_matrice{fracture_id}.under_nodes;
            ALPHA=fracture_matrice{fracture_id}.under_material;
            len=length(ALPHA);
            elem_f=[((counter+1):(counter+len))' ((counter+2):(counter+len+1))'];
            counter = counter+len+1;
        else
            elem_d=fracture_matrice{fracture_id}.above_nodes;
            ALPHA=fracture_matrice{fracture_id}.above_material;
        end
        LENGTHS = node(elem_d(:,2),:) - node(elem_d(:,1),:);
        LENGTHS = sqrt(sum(LENGTHS.^2,2));
        for j=1:2
            for i=1:2
                ii = elem_f(:,i);
                jj = elem_d(:,j);
                A = A + sparse(ii,jj,M_local(i,j)*LENGTHS.*ALPHA,N_f,N);
            end
        end
        for j=1:2
            for i=1:2
                ii = elem_d(:,i);
                jj = elem_d(:,j);
                M = M + sparse(ii,jj,M_local(i,j)*LENGTHS.*ALPHA,N,N);
            end
        end
    end
end
