%% 2 es 2
clear all
close all
clc
n = 5;
A = genera_A1(n);
% a
diag(diag(A)) % se metto diag(A,1) metto il vettore in posizione scalata
% b
Ti = tril(A) % se metto tril(A,1) prendo tutto sotto la diag +1
% triu per parte triangolare superiore
% c
detA = det(A);
rkA = rank(A);
% d
autoval = eig(A);
% e
if ne(det(A),0) % not equal
    IA = inv(A);
end
% f
autoval_inv = eig(IA);
[autoval 1./autoval_inv]
% g
% trasposta: A'
% h
N1_A=norm(A,1)
N2_A=norm(A)
max(abs(eig(A)))
% le quantità coincidono perché la matrice è simmetrica
Ninf_A=norm(A,inf)
% i
KA2 = cond(A)


