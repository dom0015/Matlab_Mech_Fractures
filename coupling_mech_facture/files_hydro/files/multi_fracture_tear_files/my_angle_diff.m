function [ c ] = my_angle_diff( a,b )
%MY_ANGLE_DIFF Summary of this function goes here
%   Detailed explanation goes here
c=a-b;
c(c<0)=c(c<0)+2*pi;

end

