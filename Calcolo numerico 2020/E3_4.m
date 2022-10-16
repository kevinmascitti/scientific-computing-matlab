%% 3 es 4
clear all
close all
clc
x0 = 0;
xN = 1;
z0 = [1;1];
f = @(x,z) [z(2); 0.1*(1-z(1)^2)*z(2)-z(1)];
N = 4;
[xE,yE] = Eulero_esplicito_sistema(f,x0,xN,N,z0);
[xH,yH] = Heun_sistema(f,x0,xN,N,z0);
plot(xE,yE,'linewidth',2);
