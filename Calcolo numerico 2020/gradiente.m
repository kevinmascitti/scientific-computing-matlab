function [x,k] = gradiente(A,b,x0,kmax,tol)
r = b-A*x0;
for k = 1:kmax
    s = A*r;
    alfa_k = r'*r/(r'*s);
    x = x0+alfa_k*r;
    r = r - alfa_k*s;
    if norm(r)/norm(b)<=tol
        break
    end
    x0 = x;
    %pause
end

