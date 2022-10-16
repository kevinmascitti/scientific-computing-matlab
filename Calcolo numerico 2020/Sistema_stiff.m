%% sistema stiff
clear all
close all
clc

A = [0 1; -100 -101];
l1 = -100; v1 = [-1;100];
l2 = -1; v2 = [1;-1];
c1 = 1/99; c2 = 100/99;
sol_esatta = @(t)[c1*exp(l1*t)*v1(1) + c2*exp(l2*t)*v2(1);...
    c1*exp(l1*t)*v1(2) + c2*exp(l2*t)*v2(2)];
x0=0;
xN=1;
xplot = linspace(x0,xN);
figure(1)
sol=sol_esatta(xplot);
plot(xplot,sol(1,:),'r')
legend('y_1')
figure(2)
plot(xplot,sol(2,:),'r')
legend('y_2')

% eulero esplicito
f = @(x,y) A*y;
y0 = [1;0];
sol=sol_esatta(xplot);

for N = [10 30 50 70 90]
    N
    [x,y] = Eulero_esplicito_sistema(f,x0,xN,N,y0);
    figure(1)
    plot(xplot,sol(1,:),'r',x,y(1,:),'--b','linewidth',2)
    legend('y_1', 'y_1 approx')
    figure(2)
    plot(xplot,sol(2,:),'r',x,y(2,:),'--b','linewidth',2)
    legend('y_2', 'y_2 approx')
    pause
end