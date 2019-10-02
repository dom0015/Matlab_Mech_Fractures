function [y] = mult_t_null(x,Z,subsizes1,subsizes2)
%BUILD_NULL Summary of this function goes here
%   Detailed explanation goes here
n=size(Z,1);
y=zeros(subsizes2(end,3),size(x,2));
for i=1:n
   y(subsizes2(i,2):subsizes2(i,3),:)=Z{i}'*x(subsizes1(i,2):subsizes1(i,3),:);
end

end