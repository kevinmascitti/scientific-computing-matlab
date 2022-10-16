function [x,k] = sezione_aurea(f,a,b,kmax,tol)

r = (sqrt(5)-1)/2;
x1 = a+(b-a)*(1-r);
x2 = a+(b-a)*r;
f1 = f(x1);
f2 = f(x2);
x = (a+b)/2;
for k = 1:kmax
    k
    hold on
    plot(x,f(x),'og','linewidth',2) 
    pause
    if f1<f2
        b = x2;
        x2 = x1;
        x1 = a +(b-a)*(1-r);
        f2 = f1;
        f1 = f(x1);
    else
        a = x1;
        x1 = x2;
        x2 = a+(b-a)*r;
        f1 = f2;
        f2 = f(x2);
    end
    x = (a+b)/2;
   if (b-a)<=tol
       break
   end
end