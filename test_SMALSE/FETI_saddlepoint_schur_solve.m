function [x,y] = FETI_saddlepoint_schur_solve(A_pinv,B_full,b,c)
%FETI_SADDLEPOINT_SCHUR_SOLVE Summary of this function goes here
%   Detailed explanation goes here


x=b;
y=c-B_full*(A_pinv*x);

x=A_pinv*x;
tmp=B_full*(A_pinv*(B_full'));

y=-tmp\y;
x=x-A_pinv*(B_full'*y);
end

