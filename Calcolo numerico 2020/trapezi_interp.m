function I = trapezi_interp(f,a,b,m)
% m = numero dei sottointervalli della partizione di [a,b]
h = (b-a)/m;
x = linspace(a,b,m+1);
y = f(x);
I = h/2*(y(1) + 2*sum(y(2:m)) + y(m+1));