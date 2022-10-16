% Equazione -u'' = f x in (a,b)
% u(a) = alfa e u(b) = beta
clear all
close all
clc

a = 0;
b = 1;
f = @(x) -20*x.^3 + 12*x.^2;
sol_esatta = @(x) x.^4.*(x-1);
alfa = 0;
beta = 0;
for k = 2:8
    % metodo diff finite con discretizzazione di [a,b]
    M = 2^k; % numero di sottointervalli: se aumenta, plot converge alla sol esatta
    h = (b-a)/M;
    x = linspace(a,b,M+1); % nodi della partizione
    % A matrice del sistema
    A = diag(2/h^2 * ones(1,M-1), 0) +...
        diag(-1/h^2 * ones(1,M-2), 1) +...
        diag(-1/h^2 * ones(1,M-2), -1);
    % termine noto con F dei nodi interni alla griglia
    F = f(x(2:M))';
    % aggiungo alfa e beta nella formula
    F(1) = F(1) + alfa/h^2;
    F(M-1) = F(M-1) +beta/h^2;
    U = A\F; % da le soluzioni dei nodi interni della partizione
    % per ottenere tutte le soluzioni aggiungo
    % le soluzioni al bordo concatenando soluzioni alfa e beta
    U = [alfa;U;beta];
    xplot = linspace(a,b,225); % serve solo un numero elevato di punti a caso
    plot(xplot,sol_esatta(xplot),'r',x,U,'o--b','linewidth',2);
    err = max(abs(sol_esatta(x) - U'))
    pause
end
% ad ogni passo l'errore viene diviso per 4 (quadratico)
