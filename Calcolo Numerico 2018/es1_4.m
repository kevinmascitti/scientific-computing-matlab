clear
close all
clc
k=1:50;
h=2.^(-k);
x=pi/4;
fe=cos(x);%derivata esatta che ci interessa
f=(sin(x+h)-sin(x))./h;
fc=2./h.*cos(x+h/2).*sin(h/2);
errrel=abs(f-fc)./abs(fc);
figure,loglog(h,errrel)
%se x è piccolo c'è cancellazione numerica, se x è grande non stiamo
%calcolando il rapporto incrementale e ci sono quindi due errori in
%competizione, risolviamo il problema usando prostaferesi