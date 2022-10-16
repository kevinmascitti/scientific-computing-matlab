%% 9 es 7
%% 1
clear all
close all
clc
% per la grafica
Pf = @(x1,x2) 100*(x2-x1.^2).^2 + (1-x1).^2;
x1 = linspace(-2,2);
x2 = x1;
x_min_ass = [1,1];
[X1,X2] = meshgrid(x1,x2);
subplot(1,2,1)
surf(X1,X2,Pf(X1,X2))
subplot(1,2,2)
contour(X1,X2,Pf(X1,X2),100)
hold on
plot(x_min_ass(1),x_min_ass(2),'*r','linewidth',2)

% f per fmincon
f = @(x) 100*(x(2)-x(1).^2).^2 +(1-x(1)).^2;
x0 = [1/2; 1/2];
A = [1 1];
b = 1;

[x,fval,exitflag,output] = fmincon(f,x0,A,b)
subplot(1,2,2)
hold on
plot(x1,1-x1,'k',x(1),x(2),'ob','linewidth',2)

% il punto di minimo vilcolato potrebbe anche coincidere cin il minimo
% assoluto,ma anche no

%% 2
clear all
close all
clc

% per la grafica
Pf = @(x1,x2) 100*(x2-x1.^2).^2 + (1-x1).^2;
x1 = linspace(-2,2);
x2 = x1;
x_min_ass = [1,1];
[X1,X2] = meshgrid(x1,x2);
subplot(1,2,1)
surf(X1,X2,Pf(X1,X2))
subplot(1,2,2)
contour(X1,X2,Pf(X1,X2),100)
hold on
plot(x_min_ass(1),x_min_ass(2),'*r','linewidth',2)

% f per fmincon
f = @(x) 100*(x(2)-x(1).^2).^2 +(1-x(1)).^2;
x0 = [1/2; 1/2];
Aeq = [1 1];
beq = 3/2;

[x,fval,exitflag,output] = fmincon(f,x0,[],[],Aeq,beq)
subplot(1,2,2)
hold on
plot(x1,3/2-x1,'k',x(1),x(2),'ob','linewidth',2)
% questo doveva stare per forza sulla retta perchè era espresso solo da una
% uguaglianza

%% 3
clear all
close all
clc
% per la grafica
Pf = @(x1,x2) 100*(x2-x1.^2).^2 + (1-x1).^2;
x1 = linspace(-2,2);
x2 = x1;
x_min_ass = [1,1];
[X1,X2] = meshgrid(x1,x2);
subplot(1,2,1)
surf(X1,X2,Pf(X1,X2))
subplot(1,2,2)
contour(X1,X2,Pf(X1,X2),100)
hold on
plot(x_min_ass(1),x_min_ass(2),'*r','linewidth',2)

% f per fmincon
f = @(x) 100*(x(2)-x(1).^2).^2 +(1-x(1)).^2;
x0 = [1/2; 1/2];
A = [1 1];
b = 1;
L = [-1 ;-inf];
U = [0; +inf];

[x,fval,exitflag,output] = fmincon(f,x0,A,b,[],[],L,U)
subplot(1,2,2)
hold on
plot(x1,1-x1,'k',x(1),x(2),'ob','linewidth',2)

%% 4
clear all
close all
clc

% per la grafica
Pf = @(x1,x2) 100*(x2-x1.^2).^2 + (1-x1).^2;
x1 = linspace(-2,2);
x2 = x1;
x_min_ass = [1,1];
[X1,X2] = meshgrid(x1,x2);
subplot(1,2,1)
surf(X1,X2,Pf(X1,X2))
subplot(1,2,2)
contour(X1,X2,Pf(X1,X2),100)
hold on
plot(x_min_ass(1),x_min_ass(2),'*r','linewidth',2)

% f per fmincon
f = @(x) 100*(x(2)-x(1).^2).^2 +(1-x(1)).^2;
x0 = [1/2; 1/2];
U = [1/2;+inf];
L = [-inf;-inf];

[x,fval,exitflag,output] = fmincon(f,x0,[],[],[],[],L,U,@vincolononlin)

% A = [1 0];
% b = 1/2;
% 
% [x,fval,exitflag,output] = fmincon(f,x0,A,b,[],[],L,U,@vincolononlin)
% riporta lo stesso risultato ma non lo stesso numero di iterazioni, perciò
% in sede d'esame non va fatto
subplot(1,2,2)
hold on

teta = linspace(0,2*pi);
xx = cos(teta);
yy = sin(teta);
plot(xx,yy,'k',x(1),x(2),'ob','linewidth',2)


%% 5
clear all
close all
clc
% per la grafica
Pf = @(x1,x2) 100*(x2-x1.^2).^2 + (1-x1).^2;
x1 = linspace(-2,2);
x2 = x1;
x_min_ass = [1,1];
[X1,X2] = meshgrid(x1,x2);
subplot(1,2,1)
surf(X1,X2,Pf(X1,X2))
subplot(1,2,2)
contour(X1,X2,Pf(X1,X2),100)
hold on
plot(x_min_ass(1),x_min_ass(2),'*r','linewidth',2)

% f per fmincon
f = @(x) 100*(x(2)-x(1).^2).^2 +(1-x(1)).^2;
x0 = [1/2; 1/2];
A = [1 1];
b = 1;
L = [-inf ;-inf];
U = [inf; +1/2];

[x,fval,exitflag,output] = fmincon(f,x0,[],[],[],[],L,U,@vincolononlin2)
subplot(1,2,2)
hold on
teta = linspace(0,2*pi);
plot(x1,1-x1.^2,'k',x(1),x(2),'ob','linewidth',2)

function  [c,ceq] = vincolononlin2(x)
    c = [];
    ceq = x(1).^2+x(2)-1;
end

function [c,ceq]= vincolononlin(x)
    c = x(1).^2+x(2).^2-1;
    ceq = [];
end