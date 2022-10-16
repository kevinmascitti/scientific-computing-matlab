%% 3 es 3
clear all
close all
clc

x0 = 0;
xN = 1;
N = 4;
y0 = [1;1];
f = @(x,y)[y(2);3*y(2)-2*y(1)];
sol_esatta = @(x) exp(x);
for N = [4 8 16 32 64];
    [x,y] = Eulero_esplicito_sistema(f,x0,xN,N,y0);
    xplot = linspace(x0,xN);
    plot(xplot,sol_esatta(xplot),'r',x,y(1,:),'--ob','linewidth',2);
    pause
end