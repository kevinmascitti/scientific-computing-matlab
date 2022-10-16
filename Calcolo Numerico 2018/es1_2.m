clear
close all
clc
n=1:16;
x=10.^(-n);%il testo diceva di interpolare in quei punti
f=(exp(x)-1)./x;
NTaylor=16;
fc=0;
for i=1:NTaylor %esegue tante iterazini quante sono le componenti di quel vettore
    %l'indice i vale per ogni iterazione il valore della componente del
    %vettore nella posizone 1,2,...ultimotermine
    fc=fc+x.^(i-1)./factorial(i);
end
errrel=abs(f-fc)./abs(fc);
figure,loglog(x,errrel)