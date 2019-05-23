function [ G ] = intersects_glue_part( intersections,materials,fracure_lengths )
%INTERSECTS_GLUE_PART Summary of this function goes here
%   Detailed explanation goes here

for i=2:length(fracure_lengths)
    intersections(intersections(:,i)>0,i)=intersections(intersections(:,i)>0,i)+sum(fracure_lengths(1:i-1));
end

n=sum(fracure_lengths)+size(intersections,1);
nn=sum(fracure_lengths);
G=sparse(n,n);
for i=1:size(intersections,1)
    temp=sparse(n,n);
    m=sum(intersections(i,:)>0);
    tf=zeros(m+1);
    tf(:,end)=-materials(i)/m;
    tf(end,:)=-materials(i)/m;
    for j=1:m
        tf(j,j)=materials(i)/m;
    end
    tf(end,end)=materials(i);
    temp([intersections(i,intersections(i,:)>0) nn+i],[intersections(i,intersections(i,:)>0) nn+i])=tf;
    G=G+temp;
end

end

