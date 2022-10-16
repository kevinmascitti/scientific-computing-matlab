%% 3 es 6
% risolvo y' = -10^3*y    x=[0,1]      y(0)=1
clear all
close all
clc
x0 = 0; 
xN = 1;
y0 = 1;
f = @(x,y) -10^3*y;
sol_esatta = @(x) x + exp(-10^3*x);
for N = [10 100 1000 10000]; % se aumenta, l'errore diminuisce
    [x,y] = Eulero_esplicito(f,x0,xN,N,y0);
    [x2,y2] = Eulero_implicito36(x0,xN,N,y0);
    % confronto grafico
    xplot = linspace(x0,xN, 201);
    subplot(1,2,1)
    semilogx(xplot, sol_esatta(xplot), 'r', x, y,...
        'ob', x2, y2, 'sg', 'linewidth', 2);
    subplot(1,2,2)
    err_rel = abs(sol_esatta(xN)-y(N+1))/abs(sol_esatta(xN));
    err_rel2 = abs(sol_esatta(xN)-y2(N+1))/abs(sol_esatta(xN));
    
    disp('   Eulero espl    Eulero impl')
    [err_rel err_rel2]
pause
end
