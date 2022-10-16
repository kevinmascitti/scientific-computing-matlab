%% 6 es 5
%% 1
clear all
close all
clc

g = @(x) -sqrt(exp(x)/2);
xplot = linspace(-1,1);
plot(xplot,g(xplot),'r',xplot,xplot,'k','linewidth',2)
x0 = -1;
kmax = 100;
tol = 1.0e-10;
[x,k] = punto_fisso(g,x0,kmax,tol)

%% 2
clear all
close all
clc

g = @(x) (2*x.^3+4*x.^2+10)./(3*x.^2+8*x);
xplot = linspace(1,2);
plot(xplot,g(xplot),'r',xplot,xplot,'k','linewidth',2)
x0 = 1.3;
kmax = 100;
tol = 1.0e-10;
[x,k] = punto_fisso(g,x0,kmax,tol)