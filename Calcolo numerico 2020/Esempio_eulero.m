%% risolviamo y'(x) = y(x)   x in [0,1]
%             y(0) = 1
clear all
close all
clc
x0 = 0;
xN = 1;
y0 = 1;
f = @(x,y) y;
sol_esatta = @(x) exp(x);
for N = [4,8,16,32,64,128,256]; % se aumenta, l'errore diminuisce
    [x,y] = Eulero_esplicito(f,x0,xN,N,y0);
    [x2,y2] = Eulero_implicito(f,x0,xN,N,y0);
    % confronto grafico
    xplot = linspace(x0,xN, 201);
    plot(xplot, sol_esatta(xplot), 'r', x, y, 'ob', x2, y2, 'sg', 'linewidth', 2);
    err_rel = abs(sol_esatta(xN)-y(N+1))/abs(sol_esatta(xN))
    err_rel2 = abs(sol_esatta(xN)-y2(N+1))/abs(sol_esatta(xN))
pause
end