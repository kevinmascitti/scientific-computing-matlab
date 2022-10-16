%% 5 es 8
clear all
close all
clc
for n = [50 100 200 400 800]
    d = ones(n,1);
    %A = spdiags([-d -d 4*d -d -d],[-n/2 -1 0 1 n/2],n,n);  % ben cond
    %A = spdiags([-d -d 3*d -d -d],[-n/2 -1 0 1 n/2],n,n);   % mal cond
    A = spdiags([-d 2*d -d],[-1 0 1],n,n);
    %full(A)
    n
    K2A = condest(A)
    %full(A)
    x_es = ones(n,1); %soluzione del sistema
    b = A*x_es;
    kmax = 200000;
    tol = 1.0e-05;
    x0 = zeros(n,1);
    % metodo del gradiente
    [x_GRAD,k_GRAD] = gradiente(A,b,x0,kmax,tol);
    err_rel_GRAD = norm(x_es-x_GRAD)/norm(x_es);
    disp(['errore realtivo GRAD = ', num2str(err_rel_GRAD), ' per numero di iterazioni = ',num2str(k_GRAD)])
    % metodo del gradiente coniugato
    [x_CG,flag,rel_res_CG,k_CG] = pcg(A,b,tol,kmax);
    err_rel_CG = norm(x_es-x_CG)/norm(x_es);
    disp(['errore realtivo PCG = ', num2str(err_rel_CG), ' per numero di iterazioni = ',num2str(k_CG)])
    % metodo del gradiente coniugato precondizionato
    L = ichol(A); % A dichiarata sparsa
    % P = L*L'  P=MN  M=L  N=L'
    [x_PCG,flag,rel_res_CG,k_PCG] = pcg(A,b,tol,kmax,L,L',x0);
    err_rel_PCG = norm(x_es-x_PCG)/norm(x_es);
    disp(['errore realtivo PCG = ', num2str(err_rel_PCG), ' per numero di iterazioni = ',num2str(k_PCG)])
    
    % con A ben cond ho circa gli stessi numeri sempre
    % con A mal cond aumentano 
end