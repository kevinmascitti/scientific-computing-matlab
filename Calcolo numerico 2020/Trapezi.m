function [x,y] = Trapezi(x0,xN,N,y0)
h = (xN-x0)/N;
x = linspace(x0,xN,N+1);
y = zeros(1,N+1);
y(1) = y0;
for n = 1:N
    y(n+1) = y(n)*(1+h/2)/(1-h/2);
end
end
