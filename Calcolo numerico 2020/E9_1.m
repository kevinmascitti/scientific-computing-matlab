% 9 es 1
% minimizzazione
% newton
clear all
close all
% a)
f = @(x) sin(x)-cos(x);
df = @(x) cos(x) + sin(x);
ddf = @(x) -sin(x) + cos(x);
fplot(f,[-2,0],'linewidth',2)
x0 = -0.5;

% b)
% f = @(x) (x-2).^4-9;
% df = @(x) 4*(x-2).^3;
% ddf = @(x) 12*(x-2).^2;
% fplot(f,[-4,8],'linewidth',2)
% x0 = 0;

% c)
% alfa = 2; x0 = 0.5;  %caso i)
%alfa = 10; x0 = 0.1; x0 = -0.1;  % caso ii)
%alfa = 10; x0 = -0.1;
%  alfa = 54; x0 = -0.1;   % caso iii)
%  alfa = 54; x0 = -0.05;
% f = @(x) x.^2 + sin(alfa*x);
% fplot(f,[-4,4],'linewidth',2)
% df = @(x) 2*x + alfa*cos(alfa*x);
% ddf = @(x) 2 - alfa^2*sin(alfa*x);
kmax = 100;
tol = 1.0e-07;
% newton
[x,k] = newton_minunc_R1(df,ddf,x0,kmax,tol,f)
hold on
plot(x,f(x),'or','linewidth',2)
minimo = f(x)