clear all
close all
clc
% matrice delle differenze finite per risolvere -u'' = f
% con u(a) = 0 e u(b) = 0
M = 10;
a = 0; 
b = 1;
h = (b-a)/M;
d = 2/h^2*ones(M-1,1);
c = -1/h^2*ones(M-2,1);
A = diag(d)+diag(c,1)+diag(c,-1);
A = sparse(A);
A = full(A);
spy(A)
numero_el_non_nulli = nnz(A)
% uso spdiags
B = [-d/2 d -d/2];
As = spdiags(B,-1:1,M-1,M-1)

whos % tipi e dimensioni delle variabili