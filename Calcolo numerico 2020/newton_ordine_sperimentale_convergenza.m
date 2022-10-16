function [x,k] = newton_ordine_sperimentale_convergenza(f,fd,x0,kmax,tol)
fr = f(x0);
for k = 1:kmax
    x = x0-f(x0)/fd(x0);
    e(k) = abs(x-x0);
    if abs(f(x))/abs(fr)<= tol
        break
    end
    x0 = x;
end

p = log(e(2:end))./log(e(1:end-1));
p'