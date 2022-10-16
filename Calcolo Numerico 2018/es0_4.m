%matrice tridiagonal ha tutti gli elementi non nulli solo sulla diagonale
%principale e quelle immeditamente vicine

%diag applicato ad una matrice tira fuori il vettore degli elementi
%diagonali, mentre se lo applico ad un vettore mi costruisce una matrice
%quadrata con tutti elementi nulli tranne la diagonale, che contiene gli
%elementi del vettore

clear
clear all
clc
%esempio
a = 1:3;
diag(a,-2);%se metto numeri diversi da 1, allora la matrice sarà più grande
%esercizio
n = 10;
x = ones(n,1);
y = ones(n-1,1); %ones crea matrici di righe e colonne dette da te

B = 5*diag(x)-diag(y,-1)+3*diag(y,1)