function [ res ] = one_fracture_nodes( nx, coo1, coo2 )
%ONE_FRACTURE Summary of this function goes here
%   Detailed explanation goes here
node1 = 1+nx*(nx+1)*coo1(1)+nx*coo1(2);
node2 = 1+nx*(nx+1)*coo2(1)+nx*coo2(2);
if coo1(1)==coo2(1)
    step = 1;
elseif coo1(2)==coo2(2)
    step = nx + 1;
else
    step = nx +2;
end
res=node1:step:node2;
end

