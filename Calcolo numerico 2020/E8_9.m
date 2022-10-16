%es_9 approssimazione di superfici
clear all
close all
clc

%R=1;
%proviamo un ellissoide
A=1;
B=4;
C=2;
theta=linspace(0,2*pi,201);
phi=linspace(0,pi,201);
[THETA,PHI]=meshgrid(theta,phi);
X=A*sin(PHI).*cos(THETA);
Y=B*sin(PHI).*sin(THETA);
Z=C*cos(PHI);
surf(X,Y,Z)

%approssimazione curva
for k=2:8
    Mt=2^k;%lo scegliamo noi
    theta_Mt=linspace(0,2*pi,Mt+1);
    Mp=2^k;%lo scegliamo noi
    phi_Mp=linspace(0,pi,Mp+1);
    [THETA_Mt,PHI_Mp]=meshgrid(theta_Mt,phi_Mp);
    
    %ora cerchiamo nodi curva
    X_Mt_Mp=A*sin(PHI_Mp).*cos(THETA_Mt);
    Y_Mt_Mp=B*sin(PHI_Mp).*sin(THETA_Mt);
    Z_Mt_Mp=C*cos(PHI_Mp);
    
        subplot(1,2,1)
    surf(X,Y,Z)
    hold on
    plot3(X_Mt_Mp,Y_Mt_Mp,Z_Mt_Mp,'or','linewidth',2)
    axis equal
    hold off
    type_approx='linear'; %'spline'
    SX=interp2(THETA_Mt,PHI_Mp,X_Mt_Mp,THETA,PHI,type_approx);
    SY=interp2(THETA_Mt,PHI_Mp,Y_Mt_Mp,THETA,PHI,type_approx);
    SZ=interp2(THETA_Mt,PHI_Mp,Z_Mt_Mp,THETA,PHI,type_approx);
        subplot(1,2,2);
    surf(SX,SY,SZ)
    axis equal
    pause
end