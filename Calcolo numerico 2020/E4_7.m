%% 4 es 7
clear all
close all
clc
a = 0;
b = 1;
E = 0.01; % oscillazioni spurie, risolvo con upwind
sol_es = @(x) (exp(x/E)-1)/(exp(1/E)-1);
xplot = linspace(a,b,401);
for k = 2:8
    M = 2^k;
    x = linspace(a,b,M+1);
    h = (b-a)/M;
    Peclet=h/(2*E)
    d = 2*E/h^2*ones(1,M-1);
    cs = (-E/h^2+1/(2*h))*ones(1,M-2);
    ci = (-E/h^2-1/(2*h))*ones(1,M-2);
    A = diag(d) + diag(cs,1) + diag(ci,-1);
    F = [zeros(M-2,1); E/h^2-1/(2*h)];
    U = A\F;
    U = [0;U;1]; % aggiungo condizioni al bordo
    plot(xplot,sol_es(xplot),'r',x,U,'--b','linewidth',2)
    err = max(abs(sol_es(x)-U'))
    % a un certo punto mi accorgo che l'errore non è più diviso per 4
    % ma quando arrivo a convergenza del metodo viene diviso per 2
    % ho approssimazione del 2 ordine per la derivata seconda
    % e del 1 ordine per la derivata prima che confonde l'ordine
    % di convergenza: lo faccio per togliere le oscillazioni spurie
    pause
end
