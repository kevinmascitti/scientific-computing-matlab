%% 4 es 3
clear all
close all
clc
a = 0;
b = 1;
y0 = exp(1) + 1;
yM = exp(1) - 1;
f = @(x) exp(x).*(1-x)-(16*pi^2+x).*sin(4*pi*x);
sol_es = @(x) exp(x) + sin(4*pi*x);
xplot = linspace(0,1);
for k=2:8
    M = 2^k;
    x = linspace(a,b,M+1);
    h = (b-a)/M;
    d = [1/h^2*ones(1,M+1) 1 1/(2*h)];
    cs = [(-2/h^2-x) 0];
    cs2 = 1/h^2*ones(1,M+1);
    A = diag(d) + diag(cs,1) + diag(cs2,2);
    A(M+2,2) = 1;
    A(M+3,1) = 1/(2*h);
    A(M+3,3) = -1/(2*h);
    A(M+3,M+1) = -1/(2*h);
    F = [f(x)'; y0; yM];
    y = A\F;
    y = y(2:M+2); % elimino nodi fittizi
    plot(xplot,sol_es(xplot),'r',x,y,'--b','linewidth',2);
    ind = [1 M/2+1 M+1];
    err = abs(sol_es(x(ind)) -y(ind)')
    pause
end