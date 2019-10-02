function [A_null]=null_elasticity_no_dirichlet(coords1,coords2)
n=2*length(coords2);

b=zeros(n,1);
b(1:2:end)=1;
b=b/norm(b);

c=zeros(n,1);
c(2:2:end)=1;
c=c/norm(c);

a=zeros(n,1);
a(1:2:end)=coords2;
a(2:2:end)=-coords1;

a=a-dot(a,b)*b;
a=a-dot(a,c)*c;
a=a/norm(a);

A_null=[a b c];

end