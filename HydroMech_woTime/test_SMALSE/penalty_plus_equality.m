%% PENALTY
addpath(genpath('files_elasticity'))
n_node=81;
[u,A_e,b_e,FREENODE,fracture_matrice,POINTS,ELEMENTS,coords1,coords2,intersections] = elasticity(n_node);
[NODE_ABOVE,NODE_UNDER] = precomp_overlap(fracture_matrice,intersections);
grad_hess=@(u)calculate_grad_hess(u, POINTS, fracture_matrice,NODE_ABOVE,NODE_UNDER);

u(FREENODE)=A_e(FREENODE,FREENODE)\b_e(FREENODE);

omega=1;
coef=1;
for j=[2 4 5]
    coef=10^j;
    for i=1:3
        fprintf('%f - %f \n',coef,omega)
        [grad_P,hess_P,~,~]=grad_hess(u);
        W=A_e+coef*hess_P;
        Aeui=A_e*u; Aeui=Aeui(FREENODE);
        A_tmp=W(FREENODE,FREENODE);
        B_tmp=(Aeui+coef*grad_P(FREENODE)-b_e(FREENODE));
        vec_tmp=A_tmp\B_tmp;
        u(FREENODE)=u(FREENODE)-omega*(vec_tmp);
    end
    figure; triplot(ELEMENTS,coords1,coords2,'g');
    hold on; triplot(ELEMENTS,coords1+u(1:2:end),coords2+u(2:2:end),'b');
    disp(norm(grad_P))
end

%% equality-constrained optimization
global switching
switching=-1;
[CROSSED_L,CROSSED_R] = find_crossed(u, POINTS, fracture_matrice,NODE_ABOVE,NODE_UNDER,0);
for j=1:1
    mat_B=calculate_B_equality(u, POINTS, fracture_matrice,NODE_ABOVE,NODE_UNDER,CROSSED_L,CROSSED_R);
    N=size(mat_B,1);
    M=[A_e(FREENODE,FREENODE) mat_B(:,FREENODE)'; mat_B(:,FREENODE) zeros(N)];
    rhs=[b_e(FREENODE); -1e-6+0*ones(N,1)];
    u_lambda=M\rhs;
    u(FREENODE)=u_lambda(1:end-N);
    figure; triplot(ELEMENTS,coords1,coords2,'g');
    hold on; triplot(ELEMENTS,coords1+u(1:2:end),coords2+u(2:2:end),'b');
    lambda=u_lambda(end-N+1:end);
    disp(lambda')
    CROSSED_L_old=CROSSED_L;
    CROSSED_R_old=CROSSED_R;
    [CROSSED_L,CROSSED_R] = find_crossed(u, POINTS, fracture_matrice,NODE_ABOVE,NODE_UNDER,1e-9);
    [CROSSED_L,CROSSED_R,changed] = combine_lambda_overlap(CROSSED_L,CROSSED_R,CROSSED_L_old,CROSSED_R_old,lambda);

%     if changed==0
%         disp(j)
%         break
%     end
    switching=-switching;
end
