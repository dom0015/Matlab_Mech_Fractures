function [Q,R] = my_QR(B)
[m,n]=size(B);
[Q,R]=qr(B');
Q=Q(:,1:m)';
R=R(1:m,:)';
end

