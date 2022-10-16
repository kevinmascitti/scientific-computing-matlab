% esercizio sistema di due equazioni non lineari
clear all
close all
clc

x = linspace(-5,5);
y = -x.^2+4*x-3;
plot(x,y,'b','linewidth',2)
x = linspace(2,5);
hold on
plot(x,sqrt(x.^2-4),'g','linewidth',2)
plot(x,-sqrt(x.^2-4),'g','linewidth',2)
x = linspace(-5,-2);
hold on
plot(x,sqrt(x.^2-4),'g','linewidth',2)
plot(x,-sqrt(x.^2-4),'g','linewidth',2)
grid on
% newton
F = @(x) [x(1).^2-x(2).^2-4; x(1).^2-4*x(1)+x(2)+3];
JF = @(x)[2*x(1), -2*x(2); 2*x(1)-4, 1];

x0 = [4;-4];
kmax = 100;
tol = 1.0e-07;
[x,k] = newton_sistemi(F,JF,x0,kmax,tol)
hold on
plot(x(1),x(2),'or','linewidth',2)
x2 = fsolve(F,x0)