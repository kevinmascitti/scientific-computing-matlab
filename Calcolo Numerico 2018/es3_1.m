clear
close all
clc

n=10;
A=hilb(n);
b=sum(A,2);
x=A\b;
xref=ones(size(x));

errrel=norm(x-xref,inf)/norm(xref,inf)
cond(A,inf)
%il numero di condizionamento � grande perch� la matrice � mal condizionata
