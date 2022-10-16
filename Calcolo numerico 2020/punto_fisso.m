function [x,k] = punto_fisso(g,x0,kmax,tol)
for k = 1:kmax
    x = g(x0);
    e(k) = abs(x-x0);
    if abs(x-x0)/abs(x) <= tol
        break
    end
    x0 = x;
end
p = log(e(2:end))./log(e(1:end-1));
p'
