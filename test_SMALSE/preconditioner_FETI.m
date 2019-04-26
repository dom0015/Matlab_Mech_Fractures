function [y] = preconditioner_FETI(x,A,A_pinv,B,r)
%PRECONDITIONER_FETI Summary of this function goes here
%   Detailed explanation goes here
m=size(A,1);
n=size(B,1);
y=zeros(m+n,1);

% y((m+1):end)=((D-1/r*W)\x((m+1):end));
% 
% y(1:m)=A_pinv*(x(1:m)-B'*y((m+1):end));

y(1:m)=x(1:m);
y((m+1):end)=x((m+1):end)-B*(A_pinv*y(1:m));

y(1:m)=A_pinv*x(1:m);

y((m+1):end)=r*y((m+1):end);

%y(1:m)=y(1:m)-A_pinv*(B'*y((m+1):end));
end

