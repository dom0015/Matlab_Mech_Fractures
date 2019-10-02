function [A_sparse] = bkdiagcell2sparse(A_cell)
%BUILD_NULL Summary of this function goes here
%   Detailed explanation goes here
n=length(A_cell);
subsizes2=zeros(n,3);
subsizes=zeros(n,3);
for i=1:n
    subsizes2(i,1)=size(A_cell{i},2);
    subsizes(i,1)=size(A_cell{i},1);
    if i==1
        subsizes2(i,2)=1;
        subsizes2(i,3)=subsizes2(i,1);
        subsizes(i,2)=1;
        subsizes(i,3)=subsizes(i,1);
    else
        subsizes2(i,2)=subsizes2(i-1,3)+1;
        subsizes2(i,3)=subsizes2(i-1,3)+subsizes2(i,1);
        subsizes(i,2)=subsizes(i-1,3)+1;
        subsizes(i,3)=subsizes(i-1,3)+subsizes(i,1);
    end
end

i_full=cell(n,1);
j_full=cell(n,1);
v_full=cell(n,1);
for i=1:n
    [ii,jj,vv]=find(A_cell{i});
    i_full{i}=ii+subsizes(i,2)-1;
    j_full{i}=jj+subsizes2(i,2)-1;
    v_full{i}=vv;
end
i_full=cat(1,i_full{:});
j_full=cat(1,j_full{:});
v_full=cat(1,v_full{:});
A_sparse=sparse(i_full,j_full,v_full,subsizes(end,3),subsizes2(end,3));
end