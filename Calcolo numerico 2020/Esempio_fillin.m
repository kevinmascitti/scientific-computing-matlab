% esempio del fill in
clear all
close all
clc
n = 100;
A = eye(n);
A(1,:) = 1;
A(:,1) = 1;
spy(A)
% fattorizzazione PA = LU
[L, U, P] = lu(A);
figure
spy(L)
title('L')
figure
spy(U)
title('U')