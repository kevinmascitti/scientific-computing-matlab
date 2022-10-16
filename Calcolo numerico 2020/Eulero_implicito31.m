function [x,y] = Eulero_implicito31(f,x0,xN,N,y0)
h = (xN-x0)/N;
x = linspace(x0,xN,N+1);
y = zeros(1,N+1);
y(1) = y0;
% y(n+1) = y(n) + h*(-y(n+1)+x(n+1)+1)
% esplicito l'incognita
% y(n+1) = (y(n)+h*x(n+1)+h) / (1+h)
for n = 1:N
    y(n+1) = (y(n)+h*x(n+1)+h) / (1+h);
end
end
