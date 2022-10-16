% 2 es 4
clear all
close all
clc
n = 4;
A = diag(2*ones(1,n), 0) + diag(-ones(1,n-1), 1) + diag(-ones(1,n-1),-1);
z = ones(n,1);
b = A*z;
[L, U, P] = lu(A);
% Ax = b -> (PA)x = Pb -> (LU)x = Pb ->
% L(Ux) = Pb -> Ux = y -> Ly = Pb
y = L\(P*b);
x = U\y;
% con comando \
x2 = A\b;
err_rel = norm(x-z)/norm(z) % migliore con PA = LU
err_rel2 = norm(x2-z)/norm(z)
KA = cond(A)
