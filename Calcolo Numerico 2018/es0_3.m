clear
close all
clc
%per adare a capo nella quadra uso ;
%per separare gli elementi della matrice uso spazi o virgole o virgole e spazi
%è uguale scrivere 1:1:6 e 1:6
a1 = 1:1:6;
a2 = 5:10;
a3 = 9:14;
a4 = 15:20;

A = [a1;a2;a3;a4]
B = A(1:end,end:-1:1) %:==1:1:end
C = A(:,2:2:end) %dalla prima riga all'ultima,prende le colonne pari
D = A(1:2:end,:) %prende solo le righe dispari
E = A([1 4 3],[5 2]) %crea una matrice con ordine ben preciso dalla matrice di A, dando le componenti
%A = A([1 4 3],[5 2]) %sovrascrive la matrice A

for ind=1:4%for si usa con notazione vettoriale secondo la variazione di ind che varia da 1 a 4, ed esegue fino a end
    d(ind)=A(ind,ind)
end %forma un vettore step by step

%diag si usa per prendere la diagonale che parte dalla posizone 1,1 di A
d = diag(A) %crea vettore colonna
