function [PRESSURE,u0,GRAD,Q,PRESSURE_diff,frac_grad,blocks] = ...
    hydro_basicSemiEvolution(D,hydro_problem,u0_old)%,divergence_diff)

const_domain = hydro_problem.const_cs_domain;
const_fracture = hydro_problem.const_cs_fracture;
const_delta_t = hydro_problem.const_delta_t;
no_fractures=hydro_problem.no_fractures;
mat_frac=hydro_problem.mat_frac;
fracture_matrice=hydro_problem.fracture_matrice;
node=hydro_problem.POINTS;
intersections=hydro_problem.intersections;
no_intersections = size(intersections,1);
alfa_inter=hydro_problem.alfa_inter;
lengths=hydro_problem.lengths;
A=hydro_problem.A;
M=hydro_problem.M;
freeNode=hydro_problem.freeNode;
b=hydro_problem.b;
u0=hydro_problem.u0;
elem=hydro_problem.ELEMENTS;

%TOCOUPLE_HANDLE Summary of this function goes here
%   Detailed explanation goes here

ALFA = cell(no_fractures,1);
MAT_FRAC = cell(no_fractures,1);
if hydro_problem.model==1
    for i=1:no_fractures
        MAT_FRAC{i} = mat_frac(i).*D{i};
        ALFA{i} = mat_frac(i)./D{i};%mat_frac(i);%./D{i};
    end
else
    for i=1:no_fractures
        MAT_FRAC{i} = mat_frac(i).*D{i}.*D{i};
        ALFA{i} = mat_frac(i).*(D{i}*0+1);
    end
end

[fracture_matrice] = a_hyd.fracture2cells_parameters( fracture_matrice,ALFA,MAT_FRAC );

N_d_node = length(A);


   
%%
[I,M]=a_hyd.interaction_matrix( fracture_matrice,node );
[F_stif,b_f,u_f,freeNode_f,F_mass] = fractures_matrix_modif( node,fracture_matrice,intersections,alfa_inter,lengths);
[MAT, freeNode, b, u0 ] = a_hyd.matrices_assembling( A, I, M, F_stif, freeNode, freeNode_f, b, b_f, u0, u_f );
% whos MAT freeNode
% disp(max(freeNode))
b = [b; 0*b_f];
u0 = [u0; 0*u_f];
MAT_TIME = [const_domain/const_delta_t*M  zeros(N_d_node,length(F_stif));
       zeros(length(F_stif),N_d_node)  const_fracture/const_delta_t*F_mass];
%MAT=MAT+MAT_TIME;

%%


MAT=MAT(freeNode,freeNode);
b0=b(freeNode);

b0=b0;%+MAT_TIME(freeNode,freeNode)*u0_old(freeNode);

x=MAT\b0;

fprintf('res_hydro=%d\n',norm(MAT*x-b0));

u0(freeNode)=x;

idx_modif = [true(length(A),1); true(length(F_stif),1)];
PRESSURE = coup.extract_pressure(u0(idx_modif),size(intersections,1),lengths);

n=length(PRESSURE);
PRESSURE_diff=cell(n,2);
for i=1:n
    tmp_f=[PRESSURE{i}(1:end-1) PRESSURE{i}(2:end)];
    PRESSURE{i}=sum(tmp_f,2)/2;
    PRESSURE_diff{i,1}=sum(tmp_f-u0(fracture_matrice{i}.above_nodes),2)/2;
    PRESSURE_diff{i,2}=sum(tmp_f-u0(fracture_matrice{i}.under_nodes),2)/2;
end

% PRESSURE_diff{1}=smoothdata(PRESSURE_diff{1},'gaussian',10);
% PRESSURE_diff{2}=smoothdata(PRESSURE_diff{2},'gaussian',10);

% figure; hold off;  plot(PRESSURE_diff{1}); hold on
% plot(PRESSURE_diff{2})

idx=length(u0)-size(intersections,1)-sum(lengths);
no_frac=length(lengths);
frac_grad=cell(no_frac,1);

for i=1:n
    %PRESSURE{i}=(temp(1:end-1)+temp(2:end))/2;
    tmp=u0((idx+1:idx+lengths(i)));
    tmp1=hydro_problem.POINTS(fracture_matrice{i}.above_nodes(:,1),:);
    tmp2=hydro_problem.POINTS(fracture_matrice{i}.above_nodes(:,2),:);
    tmp=diff(tmp);
    tmp1=tmp1-tmp2;
    tmp=tmp./sqrt(sum(tmp1.^2,2));
    frac_grad{i}=tmp;
    idx=idx+lengths(i);
end


Dirichlet_windows=hydro_problem.Dirichlet_windows;
downEdge=hydro_problem.downEdge;
rightEdge=hydro_problem.rightEdge;
upEdge=hydro_problem.upEdge;
leftEdge=hydro_problem.leftEdge;
h_elem=hydro_problem.helem;

[ Q ] = coup.extract_flow( hydro_problem.A, u0(1:size(hydro_problem.A,1)), node, elem, h_elem, Dirichlet_windows, downEdge, rightEdge, upEdge, leftEdge );

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

