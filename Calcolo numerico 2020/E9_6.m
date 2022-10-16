%% 9 es 6
% 3
% solutori matlab fminunc e fminsearch

clear all
close all
clc

f = @(x) 1/2*x(1).^2+9/2*x(2).^2;
% funzione per la grafica

fp =  @(x1,x2) 1/2*x1.^2+9/2*x2.^2;

x1 = linspace(-2,2);
x2 = x1;
[X1,X2] = meshgrid(x1,x2);
subplot(1,2,1)
surf(X1,X2,fp(X1,X2))
subplot(1,2,2)
contour(X1,X2,fp(X1,X2))

x0 = [9;1];

% solutore fminunc
[x_MU,fmu,exitflagmu,outptmu] = fminunc(f,x0)
% passando a fminunc l'espressione del gradiente di f
options = optimoptions ('fminunc','SpecifyObjectiveGradient',true);
[x_MU,fmu,exitflagmu,outptmu2] = fminunc(@f_con_grad,x0,options)
% solutore fminsearch
[x_MS,fvalms,flagms,output] = fminsearch(f,x0)

subplot(1,2,1)
hold on
plot3(x_MU(1),x_MU(2),fmu,'or','linewidth',2)

subplot(1,2,2)
hold on
plot(x_MU(1),x_MU(2),'or','linewidth',2)

function [f,Gf] = f_con_grad(x)

    f = 1/2*x(1).^2+9/2*x(2).^2;
    Gf = [x(1);9*x(2)];
end