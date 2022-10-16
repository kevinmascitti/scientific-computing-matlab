%% 4 es 8
clear all
close all
clc
a = 0;
b = 1;
E = 0.01;
sol_es = @(x) (exp(x/E)-1)/(exp(1/E)-1);
xplot = linspace(a,b,401);
for k = 2:8
    M = 2^k;
    x = linspace(a,b,M+1);
    h = (b-a)/M;
    Peclet=h/(2*E)
    d = (2*E/h^2+1/h)*ones(1,M-1);
    cs = (-E/h^2)*ones(1,M-2);
    ci = (-E/h^2-1/h)*ones(1,M-2);
    A = diag(d) + diag(cs,1) + diag(ci,-1);
    F = [zeros(M-2,1); E/h^2];
    U = A\F;
    U = [0;U;1]; % aggiungo condizioni al bordo
    plot(xplot,sol_es(xplot),'r',x,U,'--b','linewidth',2)
    err = max(abs(sol_es(x)-U'))
    pause
end
