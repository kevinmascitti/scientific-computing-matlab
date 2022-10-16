function I = simpson_interp(f,a,b,m)
% m = numero dei sottointervalli della partizione di [a,b]
h = (b-a)/(2*m);
x = linspace(a,b,2*m+1);
y = f(x);
I = h/3*(y(1) + 4*sum(y(2:2:2*m)) + 2*sum(y(3:2:2*m-1)) + y(2*m+1));