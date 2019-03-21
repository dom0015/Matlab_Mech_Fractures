function [ is_inside ] = point_inside_triangle( s, a, b, c )
%IS_POINT_INSIDE Summary of this function goes here
%   Detailed explanation goes here

as_x = s(1)-a(1);
as_y = s(2)-a(2);
s_ab = ((b(1)-a(1))*as_y-(b(2)-a(2))*as_x) > 0;
temp_ca = (c(1)-a(1))*as_y-(c(2)-a(2))*as_x > 0;

if temp_ca == s_ab
    is_inside=false;
    return;
end
temp_cb = (c(1)-b(1))*(s(2)-b(2))-(c(2)-b(2))*(s(1)-b(1)) > 0;
if temp_cb ~= s_ab
    is_inside=false;
    return;
end
    
is_inside=true;


end

