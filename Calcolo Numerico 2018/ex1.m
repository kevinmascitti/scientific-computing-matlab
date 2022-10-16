clear
close all
clc

H=hilb(10);
I=eye(10)
A=H+0.001*I;
xref=[-2:-2:-20]';
b=A*xref;
[L,U,P]=lu(A);
x=U\(L\(P*b));

Ne=norm(x-xref,inf)