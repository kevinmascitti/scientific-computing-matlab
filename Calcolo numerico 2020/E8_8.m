%esercizio 8: Approssimazione curva
clear all
close all
clc

R=1;
theta=linspace(0,2*pi,201); %discretizzazione del parametro per disegnare la curva
x=R*cos(theta);
y=R*sin(theta);


%costruire approssimante lineare a tratti e spline
for k=2:8
    M=2^k;
    theta_M=linspace(0,2*pi,M+1);
    x_M=R*cos(theta_M);
    y_M=R*sin(theta_M);
    %grafico della curva
    plot(x,y,'r','linewidth',2)
    axis equal
    hold on
    plot(x_M,y_M,'or','linewidth',2)
    x_interp=interp1(theta_M,x_M,theta);
    y_interp=interp1(theta_M,y_M,theta);
    x_spline=spline(theta_M,x_M,theta);
    y_spline=spline(theta_M,y_M,theta);
    plot(x_interp,y_interp,'--b','linewidth',2)
    plot(x_spline,y_spline,'--g','linewidth',2)
    hold off
    pause
end


