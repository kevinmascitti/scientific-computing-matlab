%esercizio1
clear all
close all
clc

%funzione runge
f=@(x) 1./(1+x.^2);
a=-5;
b=5;
for n=[5 10 15 20]                %grado del polinomio interpolante
    x=linspace(a,b,n+1);          %ascisse dati di interpolazione
    y=f(x);                       %ordinate dati di interpolazione
    c=polyfit(x,y,n);             %coeff del polinomio interpolante rispetto a base monomiale
    z=linspace(a,b,101);
    p=polyval(c,z);
    %nodi di chebichev
    xc=(b-a)/2*cos((pi*(2*[1:n+1]-1))/(2*(n+1)))+(b+a)/2 %sono le ascisse associate ai nodi di chebichev
    yc=f(xc);
    cc=polyfit(xc,yc,n);
    pc=polyval(cc,z);
    subplot(1,2,1)
    plot(z,f(z),'r',z,p,'--b',x,y,'or')
    legend('Runge','polinomio interpolante','nodi interpol')
    subplot(1,2,2)
    plot(z,f(z),'r',xc,yc,'sg',z,pc,'g')
    legend('Runge','nodi interpol','polinomio di chebichev')
    pause
end


