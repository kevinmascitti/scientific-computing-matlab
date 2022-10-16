function [x,k] = newton_minunc_R1(df,ddf,x0,kmax,tol,f)

for k = 1:kmax
    hold on
    plot(x0,f(x0),'og','linewidth',2)    
    x = x0 - df(x0)/ddf(x0);
    if abs(x-x0)/abs(x) <= tol
        break
    end
    x0 = x;  
end