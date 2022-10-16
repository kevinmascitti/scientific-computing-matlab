%% diffusione e trasporto
clear all
close all
clc
u = @(x,e) (exp(x/e)-1)/(exp(1/e)-1);
e = [10 1 0 0.1 0.01];
xplot = linspace(0,1);
for ind = 1:length(e)
    plot(xplot,u(xplot,e(ind)))
    pause
end