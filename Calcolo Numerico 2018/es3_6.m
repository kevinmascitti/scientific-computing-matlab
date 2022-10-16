clear
close all
clc

n=100;
A=rand(n);
b=sum(A,2);

tic
[Q,R]=qr(A);
XQR=R\(Q'*b);
tempoQR=toc

tic
[L,U,P]=lu(A);
XLU=U\(L\(P*b));
tempoLU=toc