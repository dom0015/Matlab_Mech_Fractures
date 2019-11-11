function [ y ] = my_shur_solver1( x,A1,A2,A3,B1,B2 )
%MY_MATMULT Summary of this function goes here
%   Detailed explanation goes here
n1=size(A1,1);
n2=size(A2,1);
n3=size(A3,1);

x1=x(1:n1);
x2=x((n1+1):(n1+n2));
x3=x((n1+n2+1):end);

y1=[x1;x2];
y2=x3;

B=[B1 B2];
A=blkdiag(A1,A2);

tmp=B'*(A3\B);
S=A-(tmp+tmp')/2;


% blok 1
y2=y2;
y1=y1-B'*(A3\y2);

%blok 2
y2=A3\y2;
y1=S\y1;

%blok 3
y1=y1;
y2=y2-A3\(B*y1);


z1=y1(1:n1);
z2=y1((n1+1):end);
y=[z1;z2;y2];

end

