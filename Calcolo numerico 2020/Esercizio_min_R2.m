% ricerca del minimo di f quadratica
clear all
close all
clc
format short e

f = @(x1,x2) 1/2*x1.^2+9/2*x2.^2;
x=linspace(-1,1);
y=linspace(-1,1);
[X,Y] = meshgrid(x,y);
subplot(1,2,1)
surf(X,Y,f(X,Y))
subplot(1,2,2)
contour(X,Y,f(X,Y),40) % fa vedere le curve di livello di f sugli assi X e Y

A = [1 0; 0 9];
b = [0;0];
kmax=100;
tol = 1.0e-06;
x0=[9;1];
[x,k]=gradiente(A,b,x0,kmax,tol);

%% esercizio fine slide
% ricerca del minimo di f quadratica
clear all
close all
clc
format short e

f = @(x1,x2) x1.^2-x1.*x2+3*x2.^2-2*x1+x2;
x=linspace(-1,2);
y=linspace(-1,1);
[X,Y] = meshgrid(x,y);
subplot(1,2,1)
surf(X,Y,f(X,Y))
subplot(1,2,2)
contour(X,Y,f(X,Y),40) % fa vedere le curve di livello di f sugli assi X e Y

A = [2 -1; -1 6];
b = [2; -1];
kmax=100;
tol = 1.0e-06;
x0=[0;0];
[x,k]=gradiente(A,b,x0,kmax,tol);
hold on
plot(x(1),x(2),'or','linewidth',2)

% newton per Rn
Gf = @(x) [2*x(1)-x(2)-2; -x(1)+6*x(2)+1];
Hf = @(x) [2 -1; -1 6];
[xN, kN] = newton_minunc_Rn(Gf,Hf,x0,kmax,tol)
subplot(1,2,2)
hold on
plot(xN(1),xN(2),'*g','linewidth',2)

%% 
clear all
clc
close all
f = @(x1,x2) x1.^4-3*x1.*x2+(x2+2).^2;
x=linspace(-4,2);
y=linspace(-8,0);
[X,Y]=meshgrid(x,y);
contour(X,Y,f(X,Y))

Gf = @(x) [4*x(1).^3-3*x(2); -3*x(1)+2*(x(2)+2)];
Hf = @(x) [12*x(1).^2 -3; -3 2];
x0 = [1;0];
kmax=100;
tol=1.0e-07;
[xN, kN] = newton_minunc_Rn(Gf,Hf,x0,kmax,tol)
hold on
plot(xN(1),xN(2),'sb','linewidth',2)