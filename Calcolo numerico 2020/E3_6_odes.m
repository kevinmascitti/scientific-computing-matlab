%% 3 es 6 con soluzione matlab
% risolvo y' = -10^3*y    x=[0,1]      y(0)=1
clear all
close all
clc
x0 = 0;
xN = 1;
y0 = 1;
f = @(x,y) -10^2*y;
sol_esatta = @(x) x + exp(-10^2*x);
xplot = linspace(x0,xN, 301);
[x45,y45] = ode45(f,[x0,xN],y0);
[x15s,y15s] = ode15s(f,[x0,xN],y0);
plot(xplot,sol_esatta(xplot),'r',x45,y45,'--b',x15s,y15s,'--g','linewidth',2)
