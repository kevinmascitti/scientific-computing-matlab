%% 4 es 1
clear all
close all
clc

a = 0;
b = pi;
f = @(x) 2*sin(x);
alfa = 1;
beta = -1;
sol_es = @(x) sin(x);
% discretizzazione
for k = 2:8
    M = 2^k;
    h = (b-a)/M;
    x = linspace(a,b,M+1);
    % matrice
    d = (2/h^2+1)*ones(1,M+1);
    c1 = [-2/h^2 -1/h^2*ones(1,M-1)];
    c2 = [-1/h^2*ones(1,M-1) -2/h^2];
    A = diag(d)+ diag(c1,1) + diag(c2,-1);
    % termine noto
    F = f(x)';
    F(1) = F(1) - 2*alfa/h;
    F(M+1) = F(M+1) +2*beta/h;
    u = A\F;
    xplot = linspace(a,b,333);
    plot(xplot,sol_es(xplot),'r',x,u,'--b','linewidth',2)
    legend('sol esatta', 'sol approx')
    err_ass(k-1) = max(abs(sol_es(x) - u')) % errore converge a 0 quadraticamente
    pause
end