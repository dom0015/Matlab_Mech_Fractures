function [ y ] = frac_sample( random_coef,grid,eig_value,eig_func )
%KL_GENERATE_RANDOM_1D Summary of this function goes here
%   Detailed explanation goes here

x=grid;
y=x*0;
for i=1:(length(random_coef))
   xi=random_coef(i);
   temp1=sqrt(eig_value(i))*eig_func{i}(x);
   temp=(temp1)*xi;
   y=y+temp;
    
end
end