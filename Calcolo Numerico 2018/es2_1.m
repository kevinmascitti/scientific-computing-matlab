clear
close all
clc
n=13;
a=-5, b=5;
f=@(dumbo)1./(1+dumbo.^2);
x=linspace(a,b,n+1);    %linspace crea n+1 nodi equispaziati da a a b
z=linspace(a,b,1001);
c=polyfit(x,f(x),n);    %calcola i coefficienti del polinomio interpolante
p=polyval(c,z);         %trova f(x) con coefficienti c e nei punti di ascissa z
plot(x,f(x),'ko',z,f(z),'b',z,p,'r');   %disegna tre funzioni, rispettivamente nera a pallini, blu e rossa