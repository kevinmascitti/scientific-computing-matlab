% interpolazione 1D
clear all
close all
clc
% funzione di runge
f = @(x) 1./(1+x.^2);
a = -5;
b = 5;
for n = 2:2:16      % grado del polinomio interpolante
    x = linspace(a,b,n+1);  % ascisse dei dati di interpolazione
    y = f(x);               % ordinate dei dati di interpolazione
    c = polyfit(x,y,n);     % coefficienti di pn rispetto alla base monomiale
    z = linspace(a,b,101);
    p = polyval(c,z);
    plot(z,f(z),'r',x,y,'or',z,p,'--b','linewidth',2)
    legend('Runge','polinomio interp')
    pause
end
