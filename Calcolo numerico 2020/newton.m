function [x,k] = newton(f,fd,x0,kmax,tol)
fr = f(x0);
for k = 1:kmax
    figure(2)
    hold on
    plot(x0,0,'og','linewidth',2)
    pause
    x = x0-f(x0)/fd(x0);
    if abs(f(x))/abs(fr)<= tol
        break
    end
    x0 = x;
end