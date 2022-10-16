clear
close all
clc

n=100;
d=ones(n,1);
dd=ones(n-1,1);

B=10*diag(d)-5*diag(dd,-1)+5*diag(dd,+1);
A=B'*B;

R=chol(A);
b=sum(A,2);

XQR=R\(R'\b);
xref=ones(size(XQR));
errrel=norm(xref-XQR)/norm(xref)

invChol=inv(R)*(inv(R))';
invref=inv(A);
norm(invref-invChol)/norm(invref)