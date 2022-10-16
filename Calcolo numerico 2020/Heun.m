function [x,y] = Heun(f, x0,xN,N,y0)
h = (xN-x0)/N;
x = linspace(x0,xN,N+1);
y = zeros(1,N+1);
y(1) = y0;
for n = 1:N
    K1 = f(x(n), y(n));
    K2 = f(x(n+1),y(n)+h*K1);
    y(n+1) = y(n)+h/2*(K1+K2);
end
end
