clear
close all
clc

x=10^(-4);
y=7-sqrt(49+x^2);
esatto=-x^2/(7+sqrt(49+x^2));
err=abs(y-esatto)/abs(esatto);

clear
close all
clc

x=10^(-12);
f1=(exp(3*x)-1)/x;
f2=3+9*x/2+27*x^2/6+81*x^3/24+243*x^4/120;
err=abs(f1-f2)/abs(f2);

x=[1,2,3,4,5;1,2,3,4,5]
sum(x)
sum(x,2)

clear
close all
clc