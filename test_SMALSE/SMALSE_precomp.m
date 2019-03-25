function [F,b0,G,u_orig,FREENODE,ELEMENTS,coords1,coords2] = SMALSE_precomp(n_node)

addpath(genpath('files_elasticity'));
[u_orig,A_e,b_e,FREENODE,fracture_matrice,POINTS,ELEMENTS,coords1,coords2,intersections] = elasticity(n_node);
[NODE_ABOVE,NODE_UNDER] = precomp_overlap(fracture_matrice,intersections);
u=u_orig;
u(FREENODE)=A_e(FREENODE,FREENODE)\b_e(FREENODE);
figure; triplot(ELEMENTS,coords1,coords2,'g');
hold on; triplot(ELEMENTS,coords1+u(1:2:end),coords2+u(2:2:end),'b');

%% equality-constrained optimization
global switching
switching=-1;
[CROSSED_L,CROSSED_R] = find_crossed(u, POINTS, fracture_matrice,NODE_ABOVE,NODE_UNDER,0);
for i=1:length(CROSSED_L)
    CROSSED_L{i}=1|CROSSED_L{i};
    CROSSED_R{i}=1|CROSSED_R{i};
end
mat_B=calculate_B_equality(u, POINTS, fracture_matrice,NODE_ABOVE,NODE_UNDER,CROSSED_L,CROSSED_R);

%% equality constrained optimization; all multiplicators "on"
N=size(mat_B,1);
M=[A_e(FREENODE,FREENODE) mat_B(:,FREENODE)'; mat_B(:,FREENODE) zeros(N)];
rhs=[b_e(FREENODE); -1e-6+0*ones(N,1)];
u_lambda=M\rhs;
u(FREENODE)=u_lambda(1:end-N);
figure; triplot(ELEMENTS,coords1,coords2,'g');
hold on; triplot(ELEMENTS,coords1+u(1:2:end),coords2+u(2:2:end),'b');

F=A_e(FREENODE,FREENODE);
b0=b_e(FREENODE,1);
G=mat_B(:,FREENODE);
