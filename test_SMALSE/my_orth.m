function [B_orht] = my_orth(B,eps_B)
%MY_ORTH Summary of this function goes here
%   Detailed explanation goes here
[m,n]=size(B);
B_orht=B;
skip_idx=true(m,1);
R=zeros(m);
for i=1:m
    B_orht(i,:)=B(i,:)/norm(B(i,:));
    for j=1:i-1
        if skip_idx(j)
            disp(dot(B_orht(i,:),B_orht(j,:)))
            B_orht(i,:)=B_orht(i,:)-dot(B_orht(i,:),B_orht(j,:))*B_orht(j,:);
        end
    end
    n_b=norm(B_orht(i,:));
    if n_b<eps_B
        fprintf('!!!!!! Not a full rank matrix! !!!!!\n');
        skip_idx(i)=false;
    else
        B_orht(i,:)=B_orht(i,:)/n_b;
    end
end
B_orht=B_orht(skip_idx,:);
end

