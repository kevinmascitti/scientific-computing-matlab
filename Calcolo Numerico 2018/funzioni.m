clear
close all
clc
%a=1;
%b=2;
%somma=f_somma(a,b)

%in un nuovo script devo scrivere SOLO
%function vegeth = f_somma(variabile 1, variabile 2)
%vegeth = variabile 1, variabile 2;
%LA FUNZIONE DEVE ESSERE SALVATA COME IL NOME DELLA FUNZIONE (es. f_somma)

x=1;
f=exp(x);
tol = 1e-10;%==10^(-10)
fe = f_TaylorExp(x,tol);
errrel = abs(f-fe)/abs(f)%errrel risulta essere minore della tolleranza
%massima data dal problema quindi ci è uscito