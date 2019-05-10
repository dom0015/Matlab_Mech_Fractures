function [PRESSURE,u0] = tocouple_handle(D,no_fractures,mat_frac,...
    fracture_matrice,node,intersections,alfa_inter,lengths,A,freeNode,b,u0)
%TOCOUPLE_HANDLE Summary of this function goes here
%   Detailed explanation goes here

ALFA = cell(no_fractures,1);
MAT_FRAC = cell(no_fractures,1);
for i=1:no_fractures
    MAT_FRAC{i} = mat_frac(i).*D{i};
    ALFA{i} = mat_frac(i)./D{i};
end

[fracture_matrice] = fracture2cells_parameters( fracture_matrice,ALFA,MAT_FRAC );

[I,M]=interaction_matrix( fracture_matrice,node );
[F,b_f,u_f,freenode_f] = fractures_matrix( node,fracture_matrice,intersections,alfa_inter,lengths);
[B, freeNode, b, u0 ] = matrices_assembling( A, I, M, F, freeNode, freenode_f, b, b_f, u0, u_f );
u0(freeNode)=B(freeNode,freeNode)\b(freeNode);

PRESSURE = extract_pressure(u0,size(intersections,1),lengths);
end

