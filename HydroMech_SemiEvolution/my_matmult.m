function [ y ] = my_matmult( x,A1,A2,A3,B1,B2 )
%MY_MATMULT Summary of this function goes here
%   Detailed explanation goes here
n1=size(A1,1);
n2=size(A2,1);
n3=size(A3,1);

x1=x(1:n1);
x2=x((n1+1):(n1+n2));
x3=x((n1+n2+1):end);

y1=A1*x1+B1'*x2;
y2=B1*x1+A2*x2+B2'*x3;
y3=B2*x2+A3*x3;
y=[y1;y2;y3];

end

