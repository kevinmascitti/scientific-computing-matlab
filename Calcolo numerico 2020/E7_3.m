% esercizio 3
clear all
close all
clc
% rappresentazione della curva
teta = linspace(0,2*pi);
r = 1;
rho = 2*r*(1+cos(teta));
x = rho.*cos(teta);
y = rho.*sin(teta);
plot(x,y,'r','linewidth',2)
axis equal
% calcolo lunghezza della curva
a = 0;
b = 2*pi;
frho = @(teta) (2*r*(1+cos(teta))).^2;
dfrho = @(teta) (-2*r*sin(teta)).^2;
f = @(teta) sqrt(frho(teta) + dfrho(teta));
I_es = 16*r;
for k = 1:8
    m(k) = 2^k;
    teta2 = linspace(0,2*pi,m(k)+1);
    rho2 = 2*r*(1+cos(teta2));
    x2 = rho2.*cos(teta2);
    y2 = rho2.*sin(teta2);
    plot(x,y,'r',x2,y2,'b','linewidth',2)
    pause
    hold off
    int_trap = trapezi(f,a,b,m(k));
    int_simp = simpson(f,a,b,m(k));
    err_trap(k) = abs(I_es - int_trap)/abs(I_es);
    err_simp(k) = abs(I_es - int_simp)/abs(I_es);    
end
[err_trap' err_simp']
% ordine sperimentale di convergenza
EOC_trap = log2(err_trap(1:end-1)./err_trap(2:end));
EOC_simp = log2(err_simp(1:end-1)./err_simp(2:end));
[EOC_trap' EOC_simp']