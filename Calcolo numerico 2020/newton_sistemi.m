function [x,k] = newton_sistemi(f,fd,x0,kmax,tol)
for k = 1:kmax
    x = x0-fd(x0)\f(x0);
    e(k) = norm(x-x0);
    if norm(x-x0)/norm(x)<= tol
        break
    end
    x0 = x;
end

p = log(e(2:end))./log(e(1:end-1));
p'
