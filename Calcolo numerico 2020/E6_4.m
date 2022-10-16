% 6 es 4
%% 1
clear all
close all
clc
 
a = 2;
f = @(x) x.^2-a;
fd = @(x) 2*x;
alfa = sqrt(a);
x0 = 0.8;
tol = 1.0e-14; 
% metodo converge velocemente anche se diminuisco la tolleranza
kmax = 100;
[x,k] = newton(f,fd,x0,kmax,tol)

%% 3 
clear all
close all
clc

f = @(x) (x-2.^(-x)).^3;
fd = @(x) 3*(x-2.^(-x)).^2.*(1 + 2.^(-x)*log(2));
xplot = linspace(-1,2);
plot(xplot,f(xplot),'r',xplot,0*xplot,'r')
x0 = 1;
tol = 1.0e-10;
kmax = 100;
[x,k] = newton(f,fd,x0,kmax,tol)

%% 4
clear all
close all
clc

f = @(x) exp(x) -2*x.^2;
fd = @(x) exp(x) -4*x;
xplot = linspace(-1,5);
plot(xplot,f(xplot),'r',xplot,0*xplot,'r')
x0 = 0.5;
tol = 1.0e-10;
kmax = 100;
[x,k] = newton(f,fd,x0,kmax,tol)
