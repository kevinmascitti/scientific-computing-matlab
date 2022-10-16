clear
close all
clc

f=@(x)1./(1+x.^2);
n=13;
x=linspace(-5,5,n+1);
z=linspace(-5,5,1001);
p=spline(x,f(x),z);

plot(x,f(x),'ko',z,f(z),'b',z,p,'r');