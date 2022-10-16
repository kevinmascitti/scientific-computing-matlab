% 2 es 3
clear all
close all
clc
n=10;
for i = 1:n
    for j = 1:n
        A(i,j) = 1/max(i,j);
    end
end
z = ones(n,1);
b = A*z;
% risolvo Ax = b con il comando \
x = A\b;
err_rel = norm(x-z)/norm(z);
% perturbazione
perturb = rand(n,1) % genera n elem casuali in 0,1
% cambiamento di variabile
% (0,1)->(a1,a2)
% t->x
% x = (a2-a1)*t + a1
a1 = -1.0e-04;
a2 = 1.0e-04;
perturb = (a2-a1)*perturb + a1;
b_perturb = b + perturb;
x_perturb = A \ b_perturb;
err_rel_x = norm(x-x_perturb) / norm(x);
err_rel_b = norm(b-b_perturb) / norm(b);
KA = cond(A)
% aumentando n, aumenta l'errore



