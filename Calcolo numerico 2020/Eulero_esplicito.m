function [x,y] = Eulero_esplicito(f,x0,xN,N,y0)
h = (xN-x0)/N;
x = linspace(x0,xN,N+1);
y = zeros(1,N+1);
y(1) = y0;
for n = 1:N
    y(n+1) = y(n) + h*f(x(n), y(n));
end
end

