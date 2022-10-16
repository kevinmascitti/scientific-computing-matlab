%es_5
clear all
close all
clc

a=-2;
b=2;
c=-2;
d=2;
f=@(x1,x2) sin(x1).*cos(x2); % exp(-10*(x1.^2+x2.^2));
nx=20;
ny=20;
x=linspace(a,b,nx+1);
y=linspace(c,d,ny+1);
z=linspace(a,b); %per la grafica
z2=linspace(c,d); %per la grafica
[Z,Z2]=meshgrid(z,z2);
%polinomio interpolante
p=0;
for i=1:nx+1
    li=polinomio_lagrange(x,i,z);
    for j=1:ny+1
        lj=polinomio_lagrange(y,j,z2);
        Lij=lj'*li;
        p=p+f(x(i),y(j))*Lij;        
    end
end
    subplot(2,2,1)
surf(Z,Z2,f(Z,Z2))
    subplot(2,2,2)
surf(Z,Z2,p)
hold on
[X,Y]=meshgrid(x,y);
plot3(X,Y,f(X,Y),'or','linewidth',2)
%aggiungo interpolante lineare a tratti
S1=interp2(X,Y,f(X,Y),Z,Z2,'linear');
    subplot(2,2,3)
surf(Z,Z2,S1)
%calcolo spline cubica
S3=interp2(X,Y,f(X,Y),Z,Z2,'spline');
    subplot(2,2,4)
surf(Z,Z2,S3)
