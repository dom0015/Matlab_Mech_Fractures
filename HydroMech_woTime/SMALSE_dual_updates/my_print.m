function [ l ] = my_print( text,l )
%PREGRESS Summary of this function goes here
%   Detailed explanation goes here
[a,b]=clear_text(l);

fprintf([a text],b{:});
l=length(text);

end

function [a,b]= clear_text(l)
a=[];
b=cell(l,1);
for i=1:l
    a=[a '%c'];
    b{i}=8;
end
end