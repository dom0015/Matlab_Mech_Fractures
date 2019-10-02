function [ frac ] = frac_from_func( node,elem,func )
%FRAC_FROM_FUNC Summary of this function goes here
%   Detailed explanation goes here
eps=0;
step=1e-5;
steps=1e4;
start=func(0);
first=node*0;
first(:,1)=node(:,1)-start(1);
first(:,2)=node(:,2)-start(2);
f_val=first(:,1).^2+first(:,2).^2;
[m,i]=min(f_val);
frac_elem=1;
frac(1)=i;
arg=0;
while arg<(1-eps)
    tmp=elem(sum(elem==frac(frac_elem),2)>0,:);
    tmp=setdiff(unique(tmp(:)),frac(frac_elem));
    tmp_coods=node(tmp,:);
    n=length(tmp);
    params=zeros(n,1);
    vals=zeros(n,1);
    for i=1:n
        [vals(i),params(i)]=find_min(tmp_coods(i,:),func,arg,step,steps);
    end
    
    [m,i]=min(vals);
    arg=params(i);
    frac_elem=frac_elem+1;
    frac(frac_elem)=tmp(i);
end
    
end

function [min_val,param]=find_min(node,func,start,step,steps)
params=start+step*(0:(steps-1));
values=sum((repmat(node,length(params),1)-func(params')).^2,2);
[m,i]=min(values);
min_val=m;
param=params(i);
end