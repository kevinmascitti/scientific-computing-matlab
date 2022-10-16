% 5 es 2
%%
clear all
clc
format short e
A = [1 -2 2; -1 1 -1; -2 -2 1];
x_es = [1;1;1];
b = A*x_es;

x0 = [0;0;0];
kmax = 400;
toll = 1.0e-7;
[xJ,kJ] = Jacobi(A,b,x0,toll,kmax);
[xG,kG] = Gauss_Seidel(A,b,x0,toll,kmax);
err_rel_J = norm(x_es-xJ)/norm(x_es);
err_rel_G = norm(x_es-xG)/norm(x_es);
disp('Jacobi')
[kJ err_rel_J]

disp('GS')
[kG err_rel_G]

%%
clear all
clc
format short e
A = [4 0 2/5; 0 5 2/5; 5/2 2 1];
x_es = [1;1;1];
b = A*x_es;

x0 = [0;0;0];
kmax = 400;
toll = 1.0e-4;
[xJ,kJ] = Jacobi(A,b,x0,toll,kmax);
[xG,kG] = Gauss_Seidel(A,b,x0,toll,kmax);
err_rel_J = norm(x_es-xJ)/norm(x_es);
err_rel_G = norm(x_es-xG)/norm(x_es);
disp('Jacobi')
[kJ err_rel_J]

disp('GS')
[kG err_rel_G]

%%
clear all
clc
format short e
A = [2 -1 1; 2 2 2; -1 -1 2];
x_es = [1;1;1];
b = A*x_es;

x0 = [0;0;0];
kmax = 400;
toll = 1.0e-4;
[xJ,kJ] = Jacobi(A,b,x0,toll,kmax);
[xG,kG] = Gauss_Seidel(A,b,x0,toll,kmax);
err_rel_J = norm(x_es-xJ)/norm(x_es);
err_rel_G = norm(x_es-xG)/norm(x_es);
disp('Jacobi')
[kJ err_rel_J]

disp('GS')
[kG err_rel_G]

%%
clear all
clc
format short e
A1 = [3 -1 0; -1 3 -1; 0 -1 3];
A2 = [0 0 -1; 0 -1 0; -1 0 0];
A = [A1 A2; A2 A1];
x_es = ones(6,1);
b = A*x_es;

x0 = zeros(6,1);
kmax = 400;
toll = 1.0e-4;
[xJ,kJ] = Jacobi(A,b,x0,toll,kmax);
[xG,kG] = Gauss_Seidel(A,b,x0,toll,kmax);
err_rel_J = norm(x_es-xJ)/norm(x_es);
err_rel_G = norm(x_es-xG)/norm(x_es);
disp('Jacobi')
[kJ err_rel_J]

disp('GS')
[kG err_rel_G]