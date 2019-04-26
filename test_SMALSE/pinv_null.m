function [A_pinv,Z,max_S] = pinv_null(A,tol)
%PINV   Pseudoinverse.
%   X = PINV(A) produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD(A) and any
%   singular values less than a tolerance are treated as zero.
%
%   PINV(A,TOL) treats all singular values of A that are less than TOL as
%   zero. By default, TOL = max(size(A)) * eps(norm(A)).
%
%   Class support for input A:
%      float: double, single
%
%   See also RANK.

%   Copyright 1984-2015 The MathWorks, Inc.

[U_orig,S_orig,V_orig] = svd(A,'econ');
s_orig = diag(S_orig);
max_S=max(s_orig);
tol=min(abs(s_orig))*tol;
r1 = sum(s_orig > tol);
V=V_orig(:,1:r1);
U=U_orig(:,1:r1);
s=s_orig(1:r1);
s = 1./s(:);
A_pinv = (V.*s.')*U';


% Orthonormal basis

Z = V_orig(:,r1+1:end);
