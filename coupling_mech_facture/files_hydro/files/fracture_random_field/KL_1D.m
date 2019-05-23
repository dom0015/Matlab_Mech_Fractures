function [ eig_value,eig_func ] = KL_1D( n,sigma,lambda )
%KL_1D Summary of this function goes here
%   Detailed explanation goes here
L=1;
[ w ] = roots_1D( n,lambda);

eig_value=zeros(n,1);
eig_func_cos=zeros(n,1);
eig_func_sin=zeros(n,1);
eig_func_x=zeros(n,1);
eig_func=cell(n,1);
for i=1:n
    eig_value(i)=(2*lambda*sigma^2)/(lambda^2*w(i)^2+1);
    eig_func_x(i)=w(i);
    eig_func_sin(i)=1/sqrt((lambda^2*w(i)^2+1)*L/2+lambda);
    eig_func_cos(i)=eig_func_sin(i)*lambda*w(i);
    eig_func{i}=@(x)eig_func_sin(i)*sin(x*eig_func_x(i))+...
        eig_func_cos(i)*cos(x*eig_func_x(i));
end


end

