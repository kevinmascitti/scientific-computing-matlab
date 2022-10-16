%% ESERCITAZIONE 1 
%% es 1a
clear
clear all
clc

x = linspace(-pi,pi,101); %vettore di punti equispaziati
f1 = @(t) sin(t);
f2 = @(t) cos(t);
plot(x,f1(x),'-r',x,f2(x),'--b','linewidth',2)
axis([-pi pi -1 1])
legend('sin(x)','cos(x)')

%% es 1b
clear all
close all
clc
x = linspace(-2,8,600);
f = @(t) exp(-t.^2).*(t>=-2 & t<0) + cos(8*t).*(t>=0 & t<2*pi) + 1*(t>=2*pi & t<=8);
plot(x,f(x),'-m','linewidth',2)
axis([-2 8 -1.2 1.2])

%% es 1c
clear all
close all
clc
x=linspace(1.0e-02,2,600);
f=@(t) t.*sin(1./t);
semilogx(x,f(x),'-c','linewidth',2)

%% es 2
clear all
close all
clc
theta=linspace(0,2*pi,500);
rho=1+0.5*cos(4*theta);
x= rho.*cos(theta);
y= rho.*sin(theta);
plot(x,y);
axis equal

%% es 3
clear all
close all
clc
t=linspace(0,10*pi,500);
a=1;
b=0.1;
x=a*cos(t);
y=a*sin(t);
z=b*t;
plot3(x,y,z);

%% es 4
clear all
close all
clc

N = 32;
x = linspace(0,1,N+1);
M = 32;
y = linspace(0,1,M+1);
[X,Y] = meshgrid(x,y);
% X ascisse di tutti i punti della griglia tensiorale
% Y ordinate di tutti i punti della matrice tensiorale
% le matrici risultanti rappresentano asc e ord della griglia
% ora valuto la funzione nei punti della griglia
f = @(x1,x2) x1.*(x1-1).*x2.*(x2-1);
f2 = @(x1,x2) x1.*(x1-1).*x2.*(x2-1).*sin(8*x1).*cos(8*x2);
mesh(X,Y,f(X,Y)) % grafico superficie
figure % apre più finestre di plot
surf(X,Y,f2(X,Y)) % grafico superficie colorato
shading flat % grafico surf senza bordo
figure
teta = linspace(0,2*pi,101);
rho = linspace(0,1,101);
[RHO, TETA] = meshgrid(rho,teta);
X = RHO.*cos(TETA);
Y = RHO.*sin(TETA);
f4 = @(x,y) exp(-10*(x.^2+y.^2));
surf(X,Y,f4(X,Y))
shading flat

%% es 8
clear all
close all
clc
format short e % forma esponenziale del risultato
f = @(x) exp(x);
df = @(x) exp(x);
x = 1;
h = 2.^(-[1:32]);
r = (f(x+h)-f(x-h))./(2*h);
err_rel = abs(df(x)-r)./abs(df(x));
err_rel' % visualizzazione in colonna
loglog(h,err_rel,'linewidth',2);
