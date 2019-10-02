function [y] = block_diag_mult(x,A,subsizes)
%BLOCK_DIAG_MULT Summary of this function goes here
%   Detailed explanation goes here
n=size(A,1);
y=sparse(size(x,1),size(x,2));
for i=1:n
    tmp1=x(subsizes(i,2):subsizes(i,3),:);
    tmp2=sparse(A{i});
    tmp=tmp2*tmp1;
    [ii,jj,vv]=find(tmp);
   %y(subsizes(i,2):subsizes(i,3),:)=tmp;
end

end

