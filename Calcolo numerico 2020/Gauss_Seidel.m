function [x,k] = Gauss_Seidel(A,b,x0,toll,kmax)

D = tril(A);
C = A-D;
% calcolo la matrice di iterazione B
B = -inv(D)*C;
rhoBG = max(abs(eigs(B)));
for k = 1:kmax
    x = D\(b-C*x0);
    err_rel = norm(x-x0)/norm(x);
    res_rel = norm(b-A*x)/norm(b);
%    if err_rel <= toll
    if res_rel <= toll
        break 
    end
    x0 = x;
end