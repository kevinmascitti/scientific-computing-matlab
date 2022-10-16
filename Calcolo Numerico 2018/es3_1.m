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
%il numero di condizionamento è grande perché la matrice è mal condizionata
