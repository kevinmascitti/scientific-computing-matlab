%esercizio 2
clear all
close all
clc

%funzione runge
f=@(x) 1./(1+x.^2);
a=-5;
b=5;
for k=1:12 %n=[5 10 15 20]                %grado del polinomio interpolante
    n=2^k;
    
    x=linspace(a,b,n+1);          %ascisse dati di interpolazione
    y=f(x);                       %ordinate dati di interpolazione
    z=linspace(a,b,2^13);
    s1=interp1(x,y,z); %spline di ordine 1 (poligonale)
    s3=spline(x,y,z);  %spline cubica
    %plot(z,f(z),'r',x,y,'or',z,s1,'--b',z,s3,':g','linewidth',2)
    err_s1(k)=norm(f(z)-s1,inf);
    err_s3(k)=norm(f(z)-s3,inf);
   % pause
end
%cerco ordine sperimentale di convergenza
EOC_s1=log2(err_s1(1:end-1)./err_s1(2:end));
EOC_s3=log2(err_s3(1:end-1)./err_s3(2:end));
[EOC_s1' EOC_s3']

