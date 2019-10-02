function [ DW ] = D_windows( i,p,n,ee )
%D_WINDOWS Summary of this function goes here
%   Detailed explanation goes here
DW = zeros(n+1,4);
DW(1,4)=p;
DW(1,2:3)=[0.0+ee     1.0-ee];
bounds = (linspace(0,1,n+1))';
for j=1:n
    DW(2:end,2)=bounds(1:end-1)+ee;
    DW(2:end,3)=bounds(2:end); DW(end,3) = 1.0-ee;
    switch i
        case 1
            DW(1,1)=3;
            DW(2:end,1)=1;
        case 2
            DW(1,1)=2;
            DW(2:end,1)=4;
        case 3
            DW(1,1)=1;
            DW(2:end,1)=3;
        case 4
            DW(1,1)=4;
            DW(2:end,1)=2;
    end

end

