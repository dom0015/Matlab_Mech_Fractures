% SMALSE David
% MPRGP

%% inputs (from elesticity with fracture)
n_node=101; % 5*k+1;
[A,b0,B,R,geometry] = SMALSE_precomp_no_dirichlet(n_node);


P_out_kerA=@(x)(x-R*(R'*x));
P_in_kerA=@(x)(R*(R'*x));
tic;
F=B*(A\P_out_kerA(B'));
toc
tmp=0*B';
tic;
for i=1:size(B,1)  
   [tmp(:,i),flag,relres,iter,resvec]=pcg(A,(P_out_kerA(B(i,:)')),1e-10,100000,diag(diag(A)));
end
toc
F=B*tmp;
imagesc(F)