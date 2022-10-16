function [x,y] = Eulero_implicito36(x0,xN,N,y0)
h = (xN-x0)/N;
x = linspace(x0,xN,N+1);
y = zeros(1,N+1);
y(1) = y0;
for n = 1:N
    y(n+1) = y(n)/(1+10^3*h);
end
end
