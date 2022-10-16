%% problema agli autovalori
clear all
close all
clc

a = 0;
b = 1;
% discretizzazione
for k = 2:8;
M = 2^k; % approssimo in maniera più accurata tutti i primi autovalori
x = linspace(a,b,M+1);
h = (b-a)/M;
% matrice A
d = 2/h^2*ones(1,M-1);
c = -1/h^2*ones(1,M-2);
A = diag(d) + diag(c,1) + diag(c,-1);
L = M-1;
Kl = pi*[1:L]'/b;
lambda_ex = Kl.^2;
% autovalori e autovettori
[U, lambda] = eig(A);
%condizioni al bordo
U = [zeros(1,M-1); U; zeros(1,M-1)];
lambda = diag(lambda);
[lambda_ex lambda]
% confronto grafico e numerico
xplot = linspace(a,b);
for ind = 1:L
    u_ex = sin(Kl(ind)*xplot);
    U_norm = sign(u_ex(2))/sign(U(2,ind))*U(:,ind)/norm(U(:,ind),inf);
    plot(xplot,u_ex,'r',x,U(:,ind),'--b','linewidth',2);
    pause
end
plot([1:L], lambda_ex, 'or', [1:L], lambda, '*b', 'linewidth', 2);
pause
end


