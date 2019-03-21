function [ POINTS, ELEMENTS, DBOUNDARY, NBOUNDARY ] = format_triangulation( p, e, t )
%FORMAT_TRIANGULATION Summary of this function goes here
%   Detailed explanation goes here
POINTS=p';
ELEMENTS=t(1:3,:)';
temp=e(1:2,e(5,:)==1);
DBOUNDARY=unique(temp(:));
NBOUNDARY=e(1:2,e(5,:)~=1)';

end

