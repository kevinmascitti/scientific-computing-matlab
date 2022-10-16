%% 5 es 4 e 6
clear all
close all
clc

n = 4000;
d = ones(n,1);
A = spdiags([-d 2*d -d],[-1 0 1],n,n);
KA = condest(A) % per matrici mal condizionate servono altri metodi
x_es = ones(n,1);
b = A*x_es;
kmax = 4000;
tol = 1.0e-5;
x0 = zeros(n,1);

[x_J,k_J] = Jacobi(A,b,x0,tol,kmax);
err_rel_J = norm(x_es-x_J)/norm(x_es);
disp(['Errore nel Jacobi = ', num2str(err_rel_J),...
    ' per numero di iterazioni = ', num2str(k_J)])

[x_GS,k_GS] = Gauss_Seidel(A,b,x0,tol,kmax);
err_rel_GS = norm(x_es-x_GS)/norm(x_es);
disp(['Errore nel Gauss-Seidel = ', num2str(err_rel_GS),...
    ' per numero di iterazioni = ', num2str(k_GS)])

[x_GRAD,k_GRAD] = gradiente(A,b,x0,kmax,tol);
err_rel_GRAD = norm(x_es-x_GRAD)/norm(x_es);
disp(['Errore nel gradiente = ', num2str(err_rel_GRAD),...
    ' per numero di iterazioni = ', num2str(k_GRAD)])