% 6 es 7
clear all
close all
clc

F = @(x) [x(1).^2+2*x(1)*x(2)+x(3);x(2).^3+x(3).^2;x(1).*x(3)-1];
JF = @(x) [2*x(1)+2*x(2), 2*x(1), 1; 0, 3*x(2)^2, 2*x(3); x(3), 0 x(1)];

x0 = [0.5;-0.5;0.1];
kmax = 20;
tol = 1.0e-10;
[x,k] = newton_sistemi(F,JF,x0,kmax,tol)

x2 = fsolve(F,x0)