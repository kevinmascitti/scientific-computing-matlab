% 6 es 1
clear all
close all
clc
f = @(x) sqrt(x.^2+1) + x.^3 + 4*x.^2 +1;
fd = @(x) x./sqrt(x.^2+1) + 3*x.^2 + 8*x;
xplot = linspace(-5,2);
figure(1)
plot(xplot,f(xplot),'r',xplot,0*xplot,'k','linewidth',2)
a = -5;
b = -2;
kmax = 100;
tol = 1.0e-07;
[xB,kB] = bisezione(f,a,b,kmax,tol)
x0 = (a+b)/2;
figure(2)
plot(xplot,f(xplot),'r',xplot,0*xplot,'k','linewidth',2)
[xN,kN] = newton(f,fd,x0,kmax,tol)
x0 = (a+b)/2;
x1 = -4;
figure(3)
plot(xplot,f(xplot),'r',xplot,0*xplot,'k','linewidth',2)
[xS,kS] = secanti(f,x0,x1,kmax,tol)