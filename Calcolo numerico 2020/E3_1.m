%% 3 es 1
% eulero implicito
clear all
close all
clc
x0 = 0;
xN = 1;
y0 = 1;
f = @(x,y) -y+x+1;
sol_esatta = @(x) x + exp(-x);
for N = [4,8,16,32,64,128,256,512]; % se aumenta, l'errore diminuisce
    [x,y] = Eulero_esplicito(f,x0,xN,N,y0);
    [x2,y2] = Eulero_implicito31(f,x0,xN,N,y0);
    [x3,y3] = Trapezi(x0,xN,N,y0);
    [x4,y4] = Heun(f,x0,xN,N,y0);
    [x5,y5] = Runge_Kutta4(f,x0,xN,N,y0);
    % confronto grafico
    xplot = linspace(x0,xN, 201);
    plot(xplot, sol_esatta(xplot), 'r', x, y, 'ob', x2, y2, 'sg',...
        x3, y3, '*c', x4, y4, 'dk',...
        x5, y5, '^m', 'linewidth', 2);
    err_rel = abs(sol_esatta(xN)-y(N+1))/abs(sol_esatta(xN));
    err_rel2 = abs(sol_esatta(xN)-y2(N+1))/abs(sol_esatta(xN));
    err_rel3 = abs(sol_esatta(xN)-y3(N+1))/abs(sol_esatta(xN));
    err_rel4 = abs(sol_esatta(xN)-y4(N+1))/abs(sol_esatta(xN));
    err_rel5 = abs(sol_esatta(xN)-y5(N+1))/abs(sol_esatta(xN));
    disp('   Eulero espl    Eulero impl       Trapz       Heun      RK4')
    [err_rel err_rel2 err_rel3 err_rel4 err_rel5]
pause
end
