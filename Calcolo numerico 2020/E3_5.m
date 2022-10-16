%% 3 es 5
%% a
clear all
close all
clc
x0 = 0;
xN = 1;
y0 = 1;
f = @(x,y) -y+x+1;
sol_esatta = @(x) x + exp(-x);
[x45,y45] = ode45(f,[x0 xN],y0);
[x23,y23] = ode23(f,[x0 xN],y0);
% vett_h = x45(2:end)-x45(1:end-1) 
% forniamo solo gli estremi  e non il tipo di intervalli
xplot = linspace(x0,xN,301);
plot(xplot,sol_esatta(xplot),'r',x45,y45,'--b',x23,y23,'--g','linewidth',2);
err_rel45 = abs(sol_esatta(xN)-y45(end))/abs(sol_esatta(xN))
err_rel23 = abs(sol_esatta(xN)-y23(end))/abs(sol_esatta(xN))

%%
clear all
close all
clc
x0 = 0;
xN = 1;
y0 = 1;
f = @(x,y) -y+x+1;
sol_esatta = @(x) x + exp(-x);
% scelgo la discretizzazione sui nodi per i metodi
N = 10;
x = linspace(x0,xN,N+1);
[x45,y45] = ode45(f,[x0 xN],y0);
[x23,y23] = ode23(f,[x0 xN],y0);
xplot = linspace(x0,xN,301);
plot(xplot,sol_esatta(xplot),'r',x45,y45,'--b',x23,y23,'--g','linewidth',2);
err_rel45 = abs(sol_esatta(xN)-y45(end))/abs(sol_esatta(xN))
err_rel23 = abs(sol_esatta(xN)-y23(end))/abs(sol_esatta(xN))
% RESULT: a parità di nodi, è migliore ode45

%%
clear all
close all
clc
x0 = 0;
xN = 1;
y0 = 1;
f = @(x,y) -y+x+1;
sol_esatta = @(x) x + exp(-x);
% definisco la tolleranza assoluta e relativa
% i solutori risolvono il problema in modo da rispettare i vincoli

options = odeset('AbsTol',1.0e-03,'RelTol',1.0e-03);
[x45,y45] = ode45(f,[x0 xN],y0,options);
[x23,y23] = ode23(f,[x0 xN],y0,options);
xplot = linspace(x0,xN,301);
plot(xplot,sol_esatta(xplot),'r',x45,y45,'--b',x23,y23,'--g','linewidth',2);
err_rel45 = abs(sol_esatta(xN)-y45(end))/abs(sol_esatta(xN))
err_rel23 = abs(sol_esatta(xN)-y23(end))/abs(sol_esatta(xN))

%% b

clear all
close all
clc

x0 = 0;
xN = 1;
y0 = [1;1];
f = @(x,y)[y(2);3*y(2)-2*y(1)];
sol_esatta = @(x) exp(x);
[x45,y45] = ode45(f,[x0 xN],y0);
[x23,y23] = ode23(f,[x0 xN],y0);
xplot = linspace(x0,xN,301);
plot(xplot,sol_esatta(xplot),'r',x45,y45(:,1),'--b','linewidth',2);
err_rel45 = abs(sol_esatta(xN)-y45(end,1))/abs(sol_esatta(xN))
err_rel23 = abs(sol_esatta(xN)-y23(end,1))/abs(sol_esatta(xN))