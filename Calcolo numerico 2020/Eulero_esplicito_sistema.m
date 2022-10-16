function [x,y] = Eulero_esplicito_sistema(f,x0,xN,N,y0)
h = (xN-x0)/N;
x = linspace(x0,xN,N+1);
m = size(y0,1); % numero di condizioni del sistema
y = zeros(m,N+1); % ora è una matrice e non un vetore riga
y(:,1) = y0;
for n = 1:N
    y(:,n+1) = y(:,n) + h*f(x(n), y(:,n));
end
end
