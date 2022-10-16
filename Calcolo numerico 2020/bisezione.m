function [x,k] = bisezione(f,a,b,kmax,tol)
fa = f(a);
% fb = f(b); % inutile
x = (a+b)/2;
% per il criterio d'arresto basato sulla tolleranza relativa sulla funzione
% f
fr = f(x);
fx = f(x);
for k = 1:kmax
    figure(1)
    hold on
    plot(x,0,'ob','linewidth',2)
    pause
    if abs(fx)/abs(fr)<=tol
        break
    else
        if fa*fx<0
            b = x;
        else
            a = x;
            fa = fx;
        end
    end
    x = (a+b)/2;
    fx = f(x);
end