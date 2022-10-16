
clear%elimina tutto ciò che è in memoria
close all%chiude tutte le immagini
clc%cancella tutto

%vettore dentro le parentesi quadre
%se uso ; esegue il comando ma non lo mostra
x = [1:-0.1:0];%parte da 1 e toglie 0.1 fino ad arrivare a 0
% x([1 4 3])
% x([1:2:7 10]) = zeros(1,5) zeros crea una matrice di 0 di m righe e n
% colonne
%x([1 2 5]) = [0.5*ones(1,2) -0.3]%ones crea una matrice di tutti 1
y = x(end:-1:1)%end legge la lunghezza del vettore, prende l'ultimo componente e vale quel numero, poi riscrive la matrice al contrario

