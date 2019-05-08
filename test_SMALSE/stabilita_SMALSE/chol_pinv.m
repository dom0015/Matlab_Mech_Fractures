function [X] = chol_pinv(A,tol)
%CHOL Summary of this function goes here
%   Detailed explanation goes here
[U,S,V] = svd(A,'econ');
s = diag(S);

r1 = sum(s > tol)+1;
V(:,r1:end) = [];
U(:,r1:end) = [];
s(r1:end) = [];
s = 1./sqrt(s(:));
X = (V.*s.')*U';

end

