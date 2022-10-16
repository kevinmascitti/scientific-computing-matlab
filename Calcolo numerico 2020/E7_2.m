% 7 es 2
clear all
close all
clc

f = @(x) 100./(x.^2).*sin(10./x);
a = 1;
b = 3;
z = linspace(a,b);
plot(z,f(z),'r','linewidth',2)
m = 16;
x = linspace(a,b,m+1);
hold on
plot(x,f(x),'b',x,0*x,'ob','linewidth',2)

% approssimazione dell'integrale
I_es = -1.42602475;
for k = 1:12
    m(k) = 2^k;
    I_trap = trapezi(f,a,b,m(k));
    I_simp = simpson(f,a,b,m(k));
    err_rel_trap(k) = abs(I_trap-I_es)/abs(I_es);
    err_rel_simp(k) = abs(I_simp-I_es)/abs(I_es);
end
% mediante function quad
[I_quad,N_quad] = quad(f,a,b,1.0e-7);
err_rel_quad = abs(I_es-I_quad)/abs(I_es)
N_quad
[(m+1)' err_rel_trap' (2*m+1)' err_rel_simp']
