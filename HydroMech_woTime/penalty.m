%% PENALTY
[u,A_e,b_e,FREENODE,fracture_matrice,POINTS,ELEMENTS,coords1,coords2,intersections] = elasticity();
A_P=@(u,NODE_above,NODE_under,crossed,overlap)penalty_matrix_exact(u, POINTS, fracture_matrice, intersections);
A_Pa=@(u,NODE_above,NODE_under,crossed,overlap)penalty_matrix_exact_approximated_intersection(u, POINTS, fracture_matrice, intersections);
A_Ps=@(u,NODE_above,NODE_under,crossed,overlap)penalty_matrix_separated(u, POINTS, fracture_matrice, intersections);
%A_P=@(u,NODE_above,NODE_under,crossed,overlap)penalty_matrix(u, POINTS, fracture_matrice, intersections);

u(FREENODE)=A_e(FREENODE,FREENODE)\b_e(FREENODE);
figure; triplot(ELEMENTS,coords1,coords2,'g');
hold on; triplot(ELEMENTS,coords1+u(1:2:end),coords2+u(2:2:end),'b');

% u=0*u;
omega=1;
coef=1;

%for j=[ 1 10 50 100 200 400 1000 2000]
for j=1:5
    coef=10^j;
    for i=1:10
        %coef=coef*1.5;
        
        fprintf('%f - %f \n',coef,omega)
%         [g_P,H_P,NODE_above,NODE_under,crossed_L,crossed_R]=A_P(u);
        [g_P,H_P,NODE_above,NODE_under,crossed_L,crossed_R,mat_B]=A_Pa(u);
        [g_P_,H_P_,~,~,~,~]=A_Ps(u);
        %[H_P,NODE_above,NODE_under,crossed_L,crossed_R]=A_Ps(u);
        %[A_Pu,NODE_above,NODE_under,crossed]=A_P(u);
        %figure; plot(sort(eig(H_P)));
        W=A_e+coef*H_P;
%         hold on; plot(sort(eig(W(FREENODE,FREENODE))))
        Aeui=A_e*u; Aeui=Aeui(FREENODE);
        %Aeui_=A_e(FREENODE,FREENODE)*u(FREENODE);
%         u(FREENODE)=u(FREENODE)+omega*(W(FREENODE,FREENODE)\b_e(FREENODE)-u(FREENODE));
        u(FREENODE)=u(FREENODE)-omega*(W(FREENODE,FREENODE)\(Aeui+coef*g_P(FREENODE)-b_e(FREENODE)));
        idx=union(NODE_above,NODE_under);
        idx=union(2*idx-1,2*idx);
%         figure; imagesc(H_P(idx,idx))
%         figure; imagesc(H_P_(idx,idx))
    end
        figure; triplot(ELEMENTS,coords1,coords2,'g');
        hold on; triplot(ELEMENTS,coords1+u(1:2:end),coords2+u(2:2:end),'b');
    res=b_e(FREENODE)-W(FREENODE,FREENODE)*u(FREENODE);
    disp(norm(res))

end

temp1=u(1:2:end);
temp2=u(2:2:end);
u_reshaped=[temp1 temp2]+POINTS;

% plot(u_reshaped(NODE_above(crossed_L),1),u_reshaped(NODE_above(crossed_L),2),'r*')
% plot(u_reshaped(NODE_under(crossed_L),1),u_reshaped(NODE_under(crossed_L),2),'m*')
% 
% plot(u_reshaped(NODE_above(crossed_R),1),u_reshaped(NODE_above(crossed_R),2),'r*')
% plot(u_reshaped(NODE_under(crossed_R),1),u_reshaped(NODE_under(crossed_R),2),'m*')

% plot(u_reshaped(NODE_above(crossed),1),u_reshaped(NODE_above(crossed),2),'r*')
% plot(u_reshaped(NODE_under(crossed),1),u_reshaped(NODE_under(crossed),2),'m*')

%% equality-constrained optimization
N=size(mat_B,1);
M=[A_e(FREENODE,FREENODE) mat_B(:,FREENODE)'; mat_B(:,FREENODE) zeros(N)];
rhs=[b_e(FREENODE); eps*ones(N,1)];
u_lambda=M\rhs;
u(FREENODE)=u_lambda(1:end-N);
lambda=u_lambda(end-N+1:end);
figure; triplot(ELEMENTS,coords1,coords2,'g');
hold on; triplot(ELEMENTS,coords1+u(1:2:end),coords2+u(2:2:end),'b');
disp(lambda)

% [CROSSED_L,CROSSED_R] = calculate_overlap(u, POINTS, fracture_matrice,intersections);