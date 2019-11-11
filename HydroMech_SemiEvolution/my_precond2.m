function [ y ] = my_precond2( x,A1,A2,A3,B1,B2 )
%MY_MATMULT Summary of this function goes here
%   Detailed explanation goes here
n1=size(A1,1);
n2=size(A3,1);
n3=size(A2,1);

x1=x(1:n1);
x3=x((n1+1):(n1+n2));
x2=x((n1+n2+1):end);

B=[B1 B2'];
A=blkdiag(A1,A2);
ia3d=spdiags(1./diag(A3),0,n2,n2);
S=A-blkdiag(B1'*inv(A3)*B1,B2*inv(A3)*B2');

y3=A3\x3;
tmp=S\([x1;x2]-B'*y3);

y1=tmp(1:n1);
y2=tmp((n1+1):end);
y=[y1;y3;y2];

end

