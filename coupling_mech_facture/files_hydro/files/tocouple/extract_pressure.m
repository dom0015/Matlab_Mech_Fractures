function PRESSURE = extract_pressure(u,no_intersections,lengths)
%EXTRACT_PRESSURE Summary of this function goes here
%   Detailed explanation goes here
idx=length(u)-no_intersections-sum(lengths);
no_frac=length(lengths);
PRESSURE=cell(no_frac,1);
for i=1:no_frac
    temp=u(idx+1:idx+lengths(i));
    PRESSURE{i}=(temp(1:end-1)+temp(2:end))/2;
    idx=idx+lengths(i);
end

