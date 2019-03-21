function [ tmp ] = my_atan( x,y )
%MY_ATAN Summary of this function goes here
%   Detailed explanation goes here

tmp=atan2(y,x);
tmp(tmp<0)=2*pi+tmp(tmp<0);
end

