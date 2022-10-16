clear
close all
clc

A=[1,2,3,4,5,6;7,8;9,10,11,12];
size(A); %dimensioni della matrice
A(1:2,4) %prendela colonna 4 dalla riga 1 alla 2 con iter 1
A(:,3) %prende tutta la colonna 3 con iter 1 da inizio a fine
A(1:2,:) %prende da riga 1 a 2 con iter 1, tutte le colonne
A(:,[2 4]) %prende tutte le righe, solo le colonne 2 e 4
A([2 3 3],:) %prende le righe 2 3 e 3 e tutte le colonne