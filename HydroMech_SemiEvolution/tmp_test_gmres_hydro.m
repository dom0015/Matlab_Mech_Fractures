
k_scale=1e-3;

A1=A(freeNode_u,freeNode_u)+const_domain/const_delta_t*M(freeNode_u,freeNode_u);
A2=F_stif+const_fracture/const_delta_t*F_mass;
A3=-k_scale^2*Au;
B1=-k_scale*B(:,freeNode_u);
B2=k_scale*G;

S=sqrt(abs([diag(A1);diag(A3);diag(A2)]));

A_fnc=@(x)my_matmult( x./S,A1,A3,A2,B1,B2)./S;


%schur=blkdiag(A1,A2)-blkdiag(B1'*diag(1./diag(A3))*B1,B2*diag(1./diag(A3))*B2');
inv_D=spdiags(1./diag(A3),0,length(A3),length(A3));
L1=chol(A1-B1'*inv_D*B1);
L2=chol(A2-B2*inv_D*B2');
P_fnc=@(x)my_precond( x./S,L1,L2,A3,B1,B2 )./S;

P_fnc2=@(x)my_precond2( x.*S,A1,A2,A3,B1,B2 ).*S;

tic;
for i=1:100
    [ y ] = my_shur_solver( b,A1,A2,A3,B1,B2 );
end
toc

