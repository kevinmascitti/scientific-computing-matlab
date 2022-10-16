clear
close all
clc

k=2;
n=2^k;
i=0:n;
x=-1+2*i/n;

f=@(x)(1-x.^2).^(5/2);
fp=@(x)(-5.*x.*(1-x.^2).^(3/2));

z=linspace(-1,1,100);
s=spline(x,f(x),z); %spline calcola e memorizza in s i valori che la spline
%cubica interpolante i dati e soddisfa le condizioni
sbarra=spline(x,[fp(x(1)),f(x),fp(x(end))],z);  %sbarra è la vincolata.
%nel primo e nell'ultimo nodo si richiede che la derivata della spline sia
%uguale a quella della funzione.
figure,plot(z,abs(f(z)-s),'b',z,abs(f(z)-sbarra),'r')
legend('not a knot','vincolata');

%la vincolata è sempre migliore se possibile