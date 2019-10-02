function [subsizes2] = null_subsizes(Z)
%BUILD_NULL Summary of this function goes here
%   Detailed explanation goes here
n=length(Z);
subsizes2=zeros(n,3);

for i=1:n
    subsizes2(i,1)=size(Z{i},2);
    if i==1
        subsizes2(i,2)=1;
        subsizes2(i,3)=subsizes2(i,1);
    else
        subsizes2(i,2)=subsizes2(i-1,3)+1;
        subsizes2(i,3)=subsizes2(i-1,3)+subsizes2(i,1);
    end
end
end