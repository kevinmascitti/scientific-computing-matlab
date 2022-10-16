% esempio modello
% ottimizzazione di newton
% minimizzazione
clear all
close all
clc
% f=@(x) sin(x)-cos(x);
% df=@(x) cos(x)+sin(x);
% ddf=@(x) -sin(x)+cos(x);

% disegno f per capire come scegliere x0
% hold on
% fplot(0,[-2,0],'linewidth',2);
% x0=-0.5;

% f=@(x) (x-2).^4-9;
% df=@(x) 4*(x-2).^3;
% ddf=@(x) 12*(x-2).^2;
% fplot(f,[-4,8]);
% x0=0;

alfa=53;
f=@(x) x.^2+sin(alfa*x);
df=@(x) 2*x+alfa*cos(alfa*x);
ddf=@(x) 2-alfa^2*sin(alfa*x);
fplot(f,[-4,4])
x0=1;


kmax=100;
tol=1e-07;
[x,k]=newton_ottim(df,ddf,x0,kmax,tol,f)
hold on
plot(x,f(x),'or','linewidth',2)
minimo=f(x)
