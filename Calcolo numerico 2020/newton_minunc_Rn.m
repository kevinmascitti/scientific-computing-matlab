function [x,k] = newton_minunc_Rn(Gf,Hf,x0,kmax,tol,f)

for k = 1:kmax   
    x = x0 - Hf(x0)\Gf(x0);
    if abs(x-x0)/abs(x) <= tol
        break
    end
    x0 = x;  
end