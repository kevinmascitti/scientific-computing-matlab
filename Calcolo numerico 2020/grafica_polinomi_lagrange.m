% rappresentazione grafica polinomi fondamentali di lagrange
clear all
close all
clc
a=0;
b=1;
c=0;
d=1;
nx=8;
ny=10;
x=linspace(a,b,nx+1);
y=linspace(c,d,ny+1);


z=linspace(a,b); %per la grafica
z2=linspace(c,d); %per la grafica

% polinomio di lagrange in 1D
% for i=1:nx+1
%     li=polinomio_lagrange(x,i,z);
%     plot(x,0*x,'or',z,li,'b');
%     pause
% end

[Z,Z2]=meshgrid(z,z2)
% polinomio di lagrange in 2D
for i=1:nx+1
    li=polinomio_lagrange(x,i,z);
    for j=1:ny+1
        lj=polinomio_lagrange(y,j,z2);
        Lij=lj'*li;
        mesh(Z,Z2,Lij);
        hold on
        plot3(x(i),y(j),1,'or','linewidth',3)
        hold off
        pause
    end
end


