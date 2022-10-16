%% 5 es 7
clear all
close all
clc
for n = [50 100 200 400 800]
    d = ones(n,1);
    %A = spdiags([-d -d 4*d -d -d],[-n/2 -1 0 1 n/2],n,n);  % ben cond
    A = spdiags([-d -d 3*d -d -d],[-n/2 -1 0 1 n/2],n,n);   % mal cond
    %full(A)
    n
    K2A = condest(A)
    x_es = ones(n,1); % soluzione del sistema
    b = A*x_es;
    kmax = 200000;
    tol = 1.0e-05;
    x0 = zeros(n,1);
    [x_GRAD,k_GRAD] = gradiente(A,b,x0,kmax,tol);
    err_rel_GRAD = norm(x_es-x_GRAD)/norm(x_es);
    disp(['errore realtivo GRAD = ', num2str(err_rel_GRAD), ' per numero di iterazioni = ',num2str(k_GRAD)])
    [x_PCG,flag,rel_res_PCG,k_PCG] = pcg(A,b,tol,kmax);
    err_rel_PCG = norm(x_es-x_PCG)/norm(x_es);
    disp(['errore realtivo PCG = ', num2str(err_rel_PCG), ' per numero di iterazioni = ',num2str(k_PCG)])
    disp('')
    
    % se metto 3*d in A, ho un errore molto grande
    % se raddoppio n, il condizionamento si quadruplica
    % lo stesso il metodo del gradiente, mentre il grad coniug raddoppia
    % gradiente coniugato è il migliore
end