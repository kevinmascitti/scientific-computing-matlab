function [x,k] = secanti(f,x0,x1,kmax,tol)
fr = f(x0);
for k = 1:kmax
    figure(3)
    hold on
    plot(x1,0,'om','linewidth',2)
    pause
    x = x1-f(x1)*(x1-x0)/(f(x1)-f(x0));
    if abs(f(x))/abs(fr)<= tol
        break
    end
    x0 = x1;
    x1 = x;
end
