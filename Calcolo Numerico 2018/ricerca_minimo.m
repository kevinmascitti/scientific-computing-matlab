clear
close all
clc

x=linspace(0.1,100,1000);
num=100*(1-0.01*x.^2).^2+0.02*x.^2;
den=(1-x.^2).^2+0.1*x.^2;
f=sqrt(num./den);

%figure,plot(x,f)
figure,loglog(x,f)%mi aiuta a vedere bene i max e i min che possono essere nascosti agli uomini
%==figure,plot(log10(x),log10(f)) anche se i valori
%sugli assi sono diversi d aquelli di loglog il disegno è uguale
%semilogx (disegna solo sull'asse x), semilogy (disegna solo sull'asse y),
%loglog scale logaritmiche da leggere nel file