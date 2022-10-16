clear
close all
clc

A=[1 2 3 4;
   -1 0 4 1;
   3 5 1 0;
   2 -1 0 1;
   1 1 -1 1;
   2 -1 0 3];
b=[1 2 3 4 5 6]';
r=rank(A);
[Q,R]=qr(A);
Rtilda=R(1:r,1:r);
c=Q'*b;
c1=c(1:r);
c2=c(r+1:end);
x=Rtilda\c1;
xref=A\b;
norm(x-xref)/norm(xref)





    