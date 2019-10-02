function [normx,normy] = fracture_boundary_values(fracture_elem_map,normx,normy,values,values2)
%FRACTURE_BOUNDARY_VALUES Summary of this function goes here
%   Detailed explanation goes here
n=length(fracture_elem_map.upi);
m=length(normx);
if nargin<5
    values2=values;
end
if ~isnumeric(values)
    
    x=linspace(0.5/n,1-0.5/n,n)';
    values=values(x);
    values2=values2(x);
end

normx((fracture_elem_map.upj-1)*m+fracture_elem_map.upi)=...
    normx((fracture_elem_map.upj-1)*m+fracture_elem_map.upi).*values;
normy((fracture_elem_map.upj-1)*m+fracture_elem_map.upi)=...
    normy((fracture_elem_map.upj-1)*m+fracture_elem_map.upi).*values;
normx((fracture_elem_map.downj-1)*m+fracture_elem_map.downi)=...
    normx((fracture_elem_map.downj-1)*m+fracture_elem_map.downi).*values2;
normy((fracture_elem_map.downj-1)*m+fracture_elem_map.downi)=...
    normy((fracture_elem_map.downj-1)*m+fracture_elem_map.downi).*values2;

end

