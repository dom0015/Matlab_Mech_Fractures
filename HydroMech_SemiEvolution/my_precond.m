function [ y ] = my_precond( x,L1,L2,A3,B1,B2 )
%MY_MATMULT Summary of this function goes here
%   Detailed explanation goes here
n1=size(L1,1);
n2=size(A3,1);
n3=size(L2,1);

x1=x(1:n1);
x3=x((n1+1):(n1+n2));
x2=x((n1+n2+1):end);

y3=A3\x3;
tmp=A3\y3;
y1=x1-B1'*tmp;
y2=x2-B2*tmp;

%schur=blkdiag(A1,A2)-blkdiag(B1'*diag(1./diag(A3))*B1,B2*diag(1./diag(A3))*B2');

z1=L1\(L1'\y1);
z2=L2\(L2'\y2);
z3=A3\(y3-B1*z1-B2'*z2);
y=[z1;z3;z2];

end

