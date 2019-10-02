function [PRESSURE,u0,GRAD,Q,PRESSURE_diff] = tocouple_handle(D,hydro_problem)


no_fractures=hydro_problem.no_fractures;
mat_frac=hydro_problem.mat_frac;
fracture_matrice=hydro_problem.fracture_matrice;
node=hydro_problem.POINTS;
intersections=hydro_problem.intersections;
alfa_inter=hydro_problem.alfa_inter;
lengths=hydro_problem.lengths;
A=hydro_problem.A;
freeNode=hydro_problem.freeNode;
b=hydro_problem.b;
u0=hydro_problem.u0;
elem=hydro_problem.ELEMENTS;



%TOCOUPLE_HANDLE Summary of this function goes here
%   Detailed explanation goes here

ALFA = cell(no_fractures,1);
MAT_FRAC = cell(no_fractures,1);
for i=1:no_fractures
    MAT_FRAC{i} = mat_frac(i).*D{i}.^2;
    ALFA{i} = 0*D{i}+1e-9;%mat_frac(i);%./D{i};
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
for i=1:1
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

%fprintf('res_hydro=%d\n',norm(A*u0(freeNode)-b0));
PRESSURE = extract_pressure(u0,size(intersections,1),lengths);

n=length(PRESSURE);
PRESSURE_diff=cell(n,2);
for i=1:n
    tmp_=[PRESSURE{i}(1:end-1) PRESSURE{i}(2:end)];
    PRESSURE{i}=sum(tmp_,2)/2;
    PRESSURE_diff{i,1}=sum(tmp_-u0(fracture_matrice{i}.above_nodes),2)/2;
    PRESSURE_diff{i,2}=sum(tmp_-u0(fracture_matrice{i}.under_nodes),2)/2;
end

Dirichlet_windows=hydro_problem.Dirichlet_windows;
downEdge=hydro_problem.downEdge;
rightEdge=hydro_problem.rightEdge;
upEdge=hydro_problem.upEdge;
leftEdge=hydro_problem.leftEdge;
h_elem=hydro_problem.helem;

[ Q ] = extract_flow( hydro_problem.A, u0(1:size(hydro_problem.A,1)), node, elem, h_elem, Dirichlet_windows, downEdge, rightEdge, upEdge, leftEdge );

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

