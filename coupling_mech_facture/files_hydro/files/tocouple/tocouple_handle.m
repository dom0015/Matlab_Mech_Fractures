function [PRESSURE,u0,GRAD] = tocouple_handle(D,no_fractures,mat_frac,...
    fracture_matrice,node,intersections,alfa_inter,lengths,A,freeNode,b,u0,elem)
%TOCOUPLE_HANDLE Summary of this function goes here
%   Detailed explanation goes here

ALFA = cell(no_fractures,1);
MAT_FRAC = cell(no_fractures,1);
for i=1:no_fractures
    MAT_FRAC{i} = mat_frac(i).*D{i}.^2;
    ALFA{i} = 0*D{i}+mat_frac(i);%./D{i};
end

[fracture_matrice] = fracture2cells_parameters( fracture_matrice,ALFA,MAT_FRAC );

[I,M]=interaction_matrix( fracture_matrice,node );
[F,b_f,u_f,freenode_f] = fractures_matrix( node,fracture_matrice,intersections,alfa_inter,lengths);
[B, freeNode, b, u0 ] = matrices_assembling( A, I, M, F, freeNode, freenode_f, b, b_f, u0, u_f );
%u0(freeNode)=B(freeNode,freeNode)\b(freeNode);


A=B(freeNode,freeNode);
b0=b(freeNode);
b=b0;
x=0*b;
y=x;
for i=1:3
b=b-A*x;
x=A\b;
y=y+x;

%fprintf('res_hydro=%d\n',norm(A*y-b0));
end
u0(freeNode)=y;
% precond=@(x)(B(freeNode,freeNode)\x)+0.01*x;
% [x,flag,relres,iter,resvec]=gmres(B(freeNode,freeNode),b(freeNode),2,1e-16,200,precond);
% figure(103);plot(resvec)
% set(gca,'YScale','log')
% u0(freeNode)=x;

fprintf('res_hydro=%d\n',norm(A*u0(freeNode)-b0));
PRESSURE = extract_pressure(u0,size(intersections,1),lengths);

coord1=node(:,1); coord2=node(:,2);
a1=coord1(elem); a2=coord2(elem);
a11=a1(:,1); a12=a1(:,2); a13=a1(:,3); 
a21=a2(:,1); a22=a2(:,2); a23=a2(:,3); 
u=u0(elem);
u1=u(:,1); u2=u(:,2); u3=u(:,3);
GRAD1=-(a21.*u2 - a22.*u1 - a21.*u3 + a23.*u1 + a22.*u3 - a23.*u2)./(a11.*a22 - a12.*a21 - a11.*a23 + a13.*a21 + a12.*a23 - a13.*a22);
GRAD2= (a11.*u2 - a12.*u1 - a11.*u3 + a13.*u1 + a12.*u3 - a13.*u2)./(a11.*a22 - a12.*a21 - a11.*a23 + a13.*a21 + a12.*a23 - a13.*a22);
GRAD=[GRAD1 GRAD2];
 
end

