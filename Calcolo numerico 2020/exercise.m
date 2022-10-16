%%
clear all
clc
close all

n=256;
alfa=2;
x=alfa*ones(1,n);
c1=ones(1,n-1);
A=diag(x)+diag(c1,1)+diag(c1,-1);
A(1,n)=0.01;
A(n,1)=0.01;

x=ones(1,n)';
b=A*x;
pcg(A,b,1.0e-06,n)

L=ichol(sparse(A));
P = L * L';
pcg(sparse(A), b, 1e-6, n, P)
%%
clear all
clc
close all

A=[-1,-4,-8;2,6,-8;-9,4,-9];
b=[-10;-2;3];
x0=zeros(1,3)';
D=diag(diag(A));
C=A-D;
B=-inv(D)*C;
rho=max(abs(eigs(B)));
for k=1:9
    x=D\(b-C*x0);
    x0=x;
end
x
%%
clear all
clc
close all

f = @(x) x^4*cos(pi*x);
df = @(x) 4*x^3*cos(pi*x) - pi*x^4*sin(pi*x);
x0 = 2;
for k=1:6
    x=x0-f(x0)/df(x0);
    x0=x;
end
x0

%%
clear all
clc
close all

f = @(x) 1/(x^2+sin(x)+1);
x0 = 0;
for k=1:100
    x=f(x0);
    e(k)=abs(x-x0);
    x0=x;
end
abs(e(10))/abs(e(9))
p=log(e(2:end))./log(e(1:end-1));
p'

%%
clear all
clc
close all

f = @(x) exp(-x^2/3);
x0=1;
for k=1:5
    x=f(x0);
    x0=x;
end
x0;

%%
clear all
clc
close all

f = @(x) log(x)-2;
df = @(x) 1/x;
x0=3;
for k=1:5
    x=x0-f(x0)/df(x0);
    x0=x;
end
x0

%%
clear all
clc
close all

f = @(x) x^3-2;
a=-1;
b=3;
fa=f(a);
fb=f(b);
x=(a+b)/2;
fx=f(x);
for k=1:3
    if fa*fx<0
        b=x;
    else
        a=x;
        fa=fx;
    end
    x=(a+b)/2;
    fx=f(x);
end
a
b

%%
clear all
clc
close all

f = @(x) x^4+4*x;
a=-1;
b=2;
x0=2;
x1=1;
for k=1:4
    x=x1-f(x1)*(x1-x0)/(f(x1)-f(x0));
    x0=x1;
    x1=x;
end
x

%%
clear all
clc
close all

f = @(x) (1/3)*x^3+3*x+2*sin(x)+pi;
a=-1;
b=0;
x=(a+b)/2;
fa=f(a);
fx=f(x);
for k=1:12
    if fa*fx<0
        b=x;
    else
        a=x;
        fa=fx;
    end
    x=(a+b)/2;
    fx=f(x);
end
x

%%
clear all
clc
close all

f = @(x) x^4+4*x;
df = @(x) 4*x^3+4;
a=-1;
b=4;
x0=2;
for k=1:5
    x=x0-f(x0)/df(x0);
    x0=x;
end
x0

%%
clear all
clc
close all

f = @(x) exp(-x+2);
a=3;
b=8;
m=98;
h=(b-a)/(2*m);
x=linspace(a,b,2*m+1);
y=f(x);
I=h/3*(y(1)+4*sum(y(2:2:2*m))+2*sum(y(3:2:2*m-1))+y(2*m+1))

%%
clear all
clc
close all

f = @(x) 1./(x.^4+1);
a=-4;
b=6;
x=linspace(a,b,7);
c=polyfit(x,f(x),6);
z=polyval(c,x);
plot(x,f(x), 'o', x, z, 'r')
hold on
fplot(f)

%%
clear all
clc
close all

f = @(x) x.*exp(x);
x=linspace(0,1,6)
c=polyfit(x,f(x),5);
z=polyval(c,x);
f(0.7)-polyval(c,0.7)
f(0.9)-polyval(c,0.9)
plot(x,z,'o',x,f(x),'r')

%%
clear all
clc
close all

x=[0.0,0.2,0.5,0.8,1.0];
f = @(x) (sqrt(x)-x.^2)./(x+1);
spline(x,f(x),0.43)

%%
clear all
clc
close all

x=[0,6,14];
y=[20,67,68];
polyval(polyfit(x,y,2),12)

%%
clear all
clc
close all
% n grado 5
% n+1 = nodi 6
f = @(x) exp(-1./x);
a=1;
b=2;
x=-cos(pi.*(2.*(0:5)+1)/12);
z=(b-a)/2.*x+(b+a)/2;
polyfit(z,f(z),5)

%%
clear all
clc
close all

x=[1,5,10];
y=[2,4,1];
spline(x,y,log(1.3))

%%
clear all
clc
close all

x=[-5,4,5,11];
y=[6,2,4,10];
y0=10;
yn=4;
spline(x,[y0 y yn],sqrt(1.8))

%%
clear all
clc
close all

f = @(x) cos(x);
x=linspace(0,pi,3);
polyval(polyfit(x,f(x),2),pi/2)

%%
clear all
clc
close all

f = @(x,y) sin(pi.*x).*sin(2*pi.*y);
n=6;
x = linspace(-1,1,n);
y=x;
[X,Y]=meshgrid(x,y);
fXY=f(X,Y);
Z_bilinear=interp2(X,Y,fXY,0.7,0.3,'linear');
Z_spline=interp2(X,Y,fXY,0.7,0.3,'spline');
abs(Z_bilinear-Z_spline)

%%
clear all
clc
close all

f=@(x,y) exp(-((sin(pi.*x)).^2-sin(pi.*y)).^2);
n=5;
x=linspace(-1,1);
y=x;
[X,Y]=meshgrid(x,y);
fXY=f(X,Y);
figure('name', 'Originale')
surf(X,Y,f(X,Y))
hold on
x=linspace(-1,1,n);
y=x;
[X,Y]=meshgrid(x,y);
fXY=f(X,Y);
figure();
Z=interp2(X,Y,fXY,X,Y,'linear');
surf(X,Y,Z)

%%
clear all
clc
close all

f = @(x,y) sin(pi.*x).*sin(2*pi.*y);
x=linspace(-1,1,6);
y=x;
[X,Y]=meshgrid(x,y);
x2=linspace(-1,1);
y2=x2;
[X2,Y2]=meshgrid(x2,y2);
Z=interp2(X,Y,f(X,Y),X2,Y2,'spline');
max(Z(:))

%%
clear all
clc
close all

f=@(x,y) exp(-sin(pi.*x).^2-sin(2*pi.*y).^2);
n=11;
x=linspace(-1,1);
y=x;
[X,Y]=meshgrid(x,y);
fXY=f(X,Y);
figure('name', 'Originale')
surf(X,Y,f(X,Y))
hold on
x=linspace(-1,1,n);
y=x;
[X,Y]=meshgrid(x,y);
fXY=f(X,Y);
figure();
Z=interp2(X,Y,fXY,X,Y,'linear');
surf(X,Y,Z)

%%
clear all
clc
close all

x0 = [1,1,1,1,1];

[x,fx,exit,output]=fmincon(@objfun, x0, [], [], [], [], [], [], @constr)

% function f = objfun(x)
%     f = abs(x(1))+abs(x(2))+abs(x(3))+abs(x(4))+abs(x(5));
% end
% 
% function [c, ceq] = constr(x)
%     c = [x(1)+x(3)-2, x(2)-x(4)-1];
%     ceq=[];
% end
%%
clear all
clc
close all

x0 = [-1,0];

[x,fx,exit,output]=fmincon(@objfun, x0, [], [], [], [], [], [], @constr)

% function f = objfun(x)
%     f = (x(1)-2).^2+x(2).^2+sin(pi*x(1));
% end
% 
% function [c, ceq] = constr(x)
%     c = x(1)+x(2).^2;
%     ceq=[];
% end
%%
clear all
clc
close all

x0 = -1;
f = @(x) exp(-(1./(1+x.^2)));
df = @(x) (2*exp(-1/(1+x.^2)).*x)./((1+x.^2).^2);
dff = @(x) (exp(-1./(1+x.^2)).*(2-6.*x.^4))./((1+x.^2).^4);
kmax=300;
for k=1:kmax
    x=x0-df(x0)/dff(x0);
    if abs(x-x0)<=1.0e-07
        break
    end
    x0=x;
end
x0
%%
clear all
clc
close all

x0 = [0,0];

[x,fx,exit,output]=fminunc(@objfun, x0)

% function f = objfun(x)
%     f = x(1).^2+(x(2)-1).^2+sin(pi.*x(1))+sin(pi.*x(2));
% end
% 
% function [c, ceq] = constr(x)
%     c = [x(1)+x(3)-2, x(2)-x(4)-1];
%     ceq=[];
% end
%%
clear all
clc
close all
f= @(x,y,z) x^2+4*y^2-x*y+5*z^2+x+y-z;
A=[2 -1 0; -1 8 0; 0 0 10];
b=[-1;-1;1];
x0=[0;1;0];
[x,k]=gradiente(A,b,x0,100,1.0e-07)
f(x(1),x(2),x(3))
%%
clear all
clc
close all

x0 = [5;2];
f = @(x,y) 10.*sin(x+y)+(x-3).^2+(y-2).^2-5;
df = @(x,y) [10.*cos(x+y)+2.*(x-3); 10.*cos(x+y)+2.*(y-2)];
dff = @(x,y) [-10.*sin(x+y)+2, -10.*sin(x+y); -10.*sin(x+y), -10.*sin(x+y)+2];
kmax=1000;
for k=1:kmax
    x=x0-dff(x0(1),x0(2))\df(x0(1),x0(2));
    if norm(x-x0)/norm(x)<=1.0e-07
        break
    end
    x0=x;
end
x0, k

%%
clear all
close all
clc

f = @(x) -sin(x)-12.*x.^3;
df = @(x) cos(x)-3.*x.^4+5;
y0 = 6;
h = 0.05;
x0 = 0;
N = 8;
xN = h*N + x0;
x = linspace(x0,xN,N+1);
y = zeros(1,N+1);
y(1) = y0;
for n=1:N
    y(n+1)=y(n)+h*f(x(n));
end
y(N+1)

%%
clear all
clc
close all

x0 = 0;
y0 = 3*pi;
xN = 6;
h = 0.01;
N = (xN-x0)/h;
f = @(x,y) y+5.*x;
x=linspace(x0,xN,N+1);
y=zeros(1,N+1);
y(1)=y0;
for n = 1:N
    y(n+1) = (y(n)+5*h.*x(n+1))/(1-h);
end
y(N)

%%
clear all
close all

f = @(x,y) [y(2); y(1)] ;
x=[0 2];
y=[sqrt(2) 0];
[m,n]=ode45(f,x,y);
n(end,1)

%%
clear all
clc
close all

f = @(x,y) -y+5.*x+2;
y0=3;
x0=0;
y=@(x)5.*x-3+6.*exp(-x);
h=0.15;
N=6;
xN=N*h+x0;
x=linspace(x0,xN,N+1);
y=zeros(1,N+1);
y(1)=y0;
for n=1:N
    K1=f(x(n),y(n));
    K2=f(x(n+1),y(n)+h*K1);
    y(n+1)=y(n)+h/2*(K1+K2);
end
y(N+1)

%%
clear all
clc
close all

f = @(x,y) -y.^2+x;
x0=0;
y0=sqrt(2);
xN=13;
[x,k]=ode45(f, [x0 xN], y0);
k(end,1)

%%
clear all
clc
close all

f = @(x,z) [z(2); 2.*z(2)-z(1)];
x0=0;
xN=1;
y0=[1;1];
N=100;
h=1/100;
x=linspace(0,1,N+1);
y=zeros(2,N+1);
y(:,1)=y0;
for n=1:N
    y(:,n+1)=y(:,n)+h*f(x(n),y(:,n));
end
err=abs(exp(1)-y(:,N+1))

%%
clear all
clc
close all

f = @(x,y) 6.*x.*y-4;
x0=0;
y0=4;
h=0.02;
N=9;
xN=h*N+x0;
x=linspace(x0,xN,N+1);
y=zeros(1,N+1);
y(1)=y0;
for n = 1:N
    y(n+1)=(y(n)-4*h)./(1-6*h.*x(n+1));
end
y(N+1)

%%
clear all
clc
close all

f = @(x,y) y.^2+x;
x0=0;
y0=pi;
xN=10;
h=0.1;
N=(xN-x0)/h;
x=linspace(x0,xN,N+1);
y=zeros(1,N+1);
y(1)=y0;
for n = 1:N
    y(n+1)=y(n)+h*f(x(n),y(n))
end
y(N+1)
    
%%
clear all
clc
close all

f=@(x,y) [y(2); 3*y(2)];
x0=0;
y0=[0;2];
xN=1;
N=100;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
y=zeros(2,N+1);
y(:,1)=y0;
for n = 1:N
    y(:,n+1)=y(:,n)+h.*f(x(n),y(:,n));
end
sol = abs(2/3*(exp(3)-1) - y(1,end))

%%
clear all
clc
close all

f = @(x,y) [y(2); y(2)-2.*y(1)-x.^3];
x0=0;
y0=[0;5];
xN=5;
h=0.1;
N=(xN-x0)/h;
x=linspace(x0,xN,N+1);
y=zeros(2,N+1);
y(:,1)=y0;
for n=1:N
    K1=f(x(n),y(:,n));
    K2=f(x(n+1),y(:,n)+h*K1);
    y(:,n+1)=y(:,n)+h/2*(K1+K2);
end
y(:,N+1)

%%
clear all
clc
close all

f=@(x) x.^3-x;
h=0.2;
xN=2.1;

df=@(x) (f(xN+h)-f(xN-h))/(2*h);
df(xN);


f=@(x) x.^2+3.*x-1;
h=0.6;
xN=1.5;

df=@(x) (f(xN+h)-f(xN))/(h);
df(xN);

f=@(x) x.^3+1;
h=0.4;
xN=2.3;

df=@(x) (f(xN)-f(xN-h))/(h);
df(xN)

%%
clear all
close all
clc

N=100;
x0=1;
xN=2;
h=1/100;
alfa=0;
beta=0;
x=linspace(x0,xN,N+1);
A = diag(1/h^2*ones(1,N-1),0) + diag(-1/(2*h^2)*ones(1,N-2),1) + diag(-1/(2*h^2)*ones(1,N-2),-1);
f=@(x) log(3.*x.^2+x);
F=f(x(1:N-1))';
U=A\F;
max(U)

%%
clear all
close all
clc

N=1000;
x0=-2;
xN=2;
h=4/1000;
x=linspace(x0,xN,N+1);
A = diag(6/h^2*ones(1,N-1),0) + diag(-3/(h^2)*ones(1,N-2),1) + diag(-3/(h^2)*ones(1,N-2),-1);
f=@(x) exp(sin(x));
F=f(x(1:N-1))';
F(1)=F(1)+3/h^2;
F(N-1)=F(N-1)+6/h^2;
U=A\F;
max(U)

%%
clear all
close all
clc

N=1000;
x0=0;
xN=pi;
h=pi/1000;
x=linspace(x0,xN,N+1);
A=diag(1/(h^2)*ones(1,N-1),0)+diag(-1/(2*h^2)*ones(1,N-2),1)+diag(-1/(2*h^2)*ones(1,N-2),-1);
A(1,1)=1/(2*h^2);
f=@(x) x-3.*cos(x);
F=f(x(1:N-1))';
F(N-1)=F(N-1)+1/(2*h^2);
U=A\F;
max(U)

%%
clear all
clc
close all

x0=0;
xN=2*pi;
N=10;
if N<11 && N>0
h=2*pi/N;
x=linspace(x0,xN,N+1);
d=-2/h^2*ones(1,N);
c=1/h^2*ones(1,N-1);
A=diag(d)+diag(c,1)+diag(c,-1);
A(N,N-1)=2/h^2;
f=@(x) -sin(x);
F=f(x(2:N+1))';
F(N)=F(N)-2/h;
U=A\F
end

%%
clc
clear all
close all

b=1;
for n=1:74
    a=b*(n^((-1)^(n+1)));
    b=a;
end
a
 
%%
clear all
clc
close all

A=[-2,7,0;7,6,3;-4,7,9];
b=[-1,-9,8];
D=diag(diag(A));
C=A-D;
B=-inv(D)*C;
max(abs(eigs(B)))

%%
clear all
clc
close all

a1=2;
b1=0.01;
n=256;
A=diag(a1*ones(1,n))+diag(ones(1,n-1),-1)+diag(ones(1,n-1),1);
A(1,n)=0.01;
A(n,1)=0.01;
x0=ones(1,n)';
b=A*x0;
L=ichol(sparse(A));
pcg(A,b,[],[],L,L');
pcg(A,b);

%%
clear all
clc
close all

A=[6,1,6;-4,-2,-3;4,-9,2];
b=[5;-8;-8];
x0=[0,0,0]';
D=diag(diag(A));
C=A-D;
for k=1:8
    x=D\(b-C*x0);
    x0=x;
end
x0

%%
clear all
clc
close all

n=256;
alfa=2;
d=alfa*ones(1,n);
f=-1*ones(1,n-1);
A=diag(d)+diag(f,1)+diag(f,-1);
x0=ones(1,n)';
b=A*x0;
pcg(A,b,1.0e-06,n)

%%
clear all
close all
clc

f=@(x) x^3+3*x+sin(pi*x);
df=@(x) 3*x^2+3+pi*cos(pi*x);
a=-2;
b=3;
x0=2;
kmax=5;
for k=1:kmax
    x=x0-f(x0)/df(x0);
    x0=x;
end
x0

%%
clear all
close all
clc

f=@(x) x^3*sin(x);
df=@(x) 3*x^2*sin(x)+x^3*cos(x);
x0=2;
kmax=4;
for k=1:kmax
    x=x0-f(x0)/df(x0);
    x0=x;
end
x0

%%
clear all
close all
clc

f=@(x) exp(-x^2/3);
x0=1;
for k=1:5
    x=f(x0);
    x0=x;
end
x0

%%
clear all
clc
close all
f=@(x) x^4+2*x^3-4*x-8;
a=0;
b=4;
x=(a+b)/2;
fa=f(a);
fx=f(x);
for k=1:6
    if fa*fx<0
        b=x;
    else
        a=x;
    end
    x=(a+b)/2;
    fa=f(a);
    fx=f(x);
end
x

%%
clear all
clc
close all

f=@(x) 1/(x^2+sin(x)+1);
x0=0;
for k=1:20
    x=f(x0);
    e(k)=abs(x-x0);
    x0=x;
end
e(10)/e(9)
p=log(e(2:end))./log(e(1:end-1));
p'

%%
clear all
close all
clc

f=@(x) x^3+sin(x);
x0=-2;
x1=1;
N=3;
for k=1:N
    x=x1-f(x1)*(x1-x0)/(f(x1)-f(x0));
    x0=x1;
    x1=x;
end
x1

%%
clear all
clc
close all

f=@(x) log(x)-1;
df=@(x) 1/x;
x0=3;
for k=1:4
    x=x0-f(x0)/df(x0);
    x0=x;
end
x0
    
%%
clear all
close all
clc

f=@(x) (x(1)-2).^2+x(2).^2+sin(pi.*x(1));
x0=[0;0];
fminunc(f,x0)

%%
clear all
clc
close all

f=@(x) abs(x(1)-1)+abs(x(2)+2)-1;
g=@(x,y) abs(x-1)+abs(y+2)-1;
x0=[1;1];
[x,y,exit,output] = fmincon(f,x0,[-1 1],-2);
output
x
y
[X,Y]=meshgrid(linspace(-3,3),linspace(-3,3));
surf(X,Y,g(X,Y))
hold on
plot3(1,-2,g(1,-2),'or')

%%
clear all
clc
close all

f=@(x) x(1).^2+(x(2)-1).^2+sin(pi.*x(1))+sin(pi.*x(2));
x0=[4,4];
[x,y,exit,output] = fmincon(f,x0,[],[],[],[],[0,-Inf],[Inf,0],@costraint);
x,y,output
% function [c,ceq] = costraint(x)
%     c = x(1).*x(2)+1;
%     ceq= [] ;
% end

%%
clear all
clc
close all

f=@(x) exp(-1/(1+x.^2));
x0=-1;
[x,y,exit,output]=fminunc(f,x0);
x,y,output

%%
clear all
clc
close all

f=@(x) 10.*sin(x(1)+x(2))+(x(1)-3).^2+(x(2)-2).^2-5;
x0=[0;0];
[x,y,exit,output]=fminunc(f,x0);
x,y,output
g=@(x,y) 10.*sin(x+y)+(x-3).^2+(y-2).^2-5;
[X,Y]=meshgrid(linspace(-3,3),linspace(-3,3));
surf(X,Y,g(X,Y))
hold on
plot3(0.0323,-0.9677,g(0.0323,-0.9677),'or')

%%
clear all
clc
close all

f=@(x) x(1).^2+4.*x(2).^2-x(1).*x(2)+5.*x(3).^2+x(1)+x(2)-x(3);
x0=[0;1;0];
[x,y,exit,output] = fminunc(f,x0);
x,y,output

%%
clear all
clc
close all

x=[-1,1,7,9,19];
y=[4,3,10,10,9];
spline(x,y,log(0.9))

%%
clear all
clc
close all

f=@(x) x.^2+sin(x.^3);
N=5;
x=linspace(0,1,N);
c=polyfit(x,f(x),N-1);
a=polyval(c,0.3);
b=polyval(c,0.5);

abs(f(0.3)-a)
abs(f(0.5)-b)


%%
clear all
clc
close all
    
N=6;
x=linspace(0,1,N);
f=@(x) (sin(x)-(x+1).^2)./(sqrt(x)-2);
c=spline(x,f(x),0.6)


%%
clear all
clc
close all
    
x=[4,6,8,9];
y=[6,4,5,9];
spline(x,[8 y 5], exp(0.4))

%%
clear all
close all
clc

t=[0,10,14,22,28];
v=[42,43,38,61,50];
c=polyfit(t,v,4);
z=polyval(c,24)

%%
clear all
clc
close all
    
f=@(x) cos(x);
N=3;
x=linspace(0,2*pi,N);
c=polyfit(x,f(x),N-1);
z=polyval(c,pi/2)

%%
clear all
clc
close all
    
f=@(x) sin(x);
N=5;
x=linspace(0,pi,N);
c=polyfit(x,f(x),N-1);
z=polyval(c,pi/8)

%%
clear all
clc
close all

x=[2,9,15,20];
y=[8,10,5,9];
c=polyfit(x,y,3);
z=polyval(c,exp(0.5))

%%
clear all
clc
close all

f=@(x,y) exp(-(sin(2*pi.*x)).^2-(sin(2*pi.*y)).^2);
x=linspace(-1,1);
[X,Y]=meshgrid(x,x);
figure(1)
surf(X,Y,f(X,Y))
n=7;
x=linspace(-1,1,n);
[X2,Y2]=meshgrid(x,x);
Z=interp2(X2,Y2,f(X2,Y2),X,Y,'bilinear');
figure(2)
surf(X,Y,Z)

%%
clear all
clc
close all

f=@(x,y) exp(-(sin(pi.*x)).^2-(sin(pi.*y)).^2);
x=linspace(-1,1);
[X,Y]=meshgrid(x,x);
figure(1)
surf(X,Y,f(X,Y))
n=21;
x=linspace(-1,1,n);
[X2,Y2]=meshgrid(x,x);
Z=interp2(X2,Y2,f(X2,Y2),X,Y,'bilinear');
figure(2)
surf(X,Y,Z)

%%
clear all
clc
close all

f=@(x,y) sin(pi.*x).*sin(2*pi.*y);
x=linspace(-1,1);
[X,Y]=meshgrid(x,x);
figure(1)
surf(X,Y,f(X,Y))
n=6;
x=linspace(-1,1,n);
[X2,Y2]=meshgrid(x,x);
Z=interp2(X2,Y2,f(X2,Y2),X,Y,'spline');
max(Z(:))

%%
clear all
clc
close all

f=@(x,y) sin(pi.*x).*sin(2*pi.*y);
x=linspace(-1,1);
[X,Y]=meshgrid(x,x);
% figure(1)
% surf(X,Y,f(X,Y))
n=6;
x=linspace(-1,1,n);
[X2,Y2]=meshgrid(x,x);
Z=interp2(X2,Y2,f(X2,Y2),0.7,0.3,'bilinear');
Z2=interp2(X2,Y2,f(X2,Y2),0.7,0.3,'spline');
abs(Z-Z2)

%%
clear all
clc
close all

N=12;
x0=-8;
xN=5;
h=(xN-x0)/(2*N);
x=linspace(x0,xN,2*N+1);
f=@(x) cos(3*x+3)-4*x;
y=f(x);
I=h/3*(y(1)+4*sum(y(2:2:2*N))+2*sum(y(3:2:2*N-1))+y(2*N+1))

%%
clear all
clc
close all

f=@(x) exp(-x.^2);
x0=-1;
xN=1;
N=12;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
y=f(x);
I=h/2*(y(1)+2*sum(y(2:N))+y(N+1))

%%
clear all
clc
close all

f=@(x) 3*exp(-x);
x0=2;
xN=4;
N=50;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
y=f(x);
I=h/2*(y(1)+2*sum(y(2:N))+y(N+1))

%%
clear all
close all
clc

f=@(x,y) cos(x.^2)./y.^2;
x=linspace(-pi,pi,6);
y=linspace(1,3,10);
[X,Y]=meshgrid(x,y);
[Xtot,Ytot]=meshgrid(linspace(-pi,pi),linspace(1,3));
Z1=interp2(X,Y,f(X,Y),Xtot,Ytot,'spline');
Z2=interp2(X,Y,f(X,Y),Xtot,Ytot,'linear');
max(max(abs(Z1-Z2)))

%%
clear all
clc
close all

h=0.05;
f=@(x)x.^2.*exp(-x.^2+1);
df=@(x)2.*x.*exp(-x.^2+1)-2.*x.^3.*exp(-x.^2+1);
((f(2+h)-f(2))/h-df(2))/(df(2))

%%
clear all
clc
close all

N=16;
xN=2;
x0=0;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
d=(2/h^2)*ones(1,N-1);
c1=(-1/h^2-10/h)*ones(1,N-2);
c2=(-1/h^2+10/h)*ones(1,N-2);
A=diag(d)+diag(c1,-1)+diag(c2,1);
u=@(x) 2./(exp(40)-1).*(exp(20.*x)-1);
F=zeros(1,N-1)';
F(N-1)=2/h^2-20/h;
U=A\F;
U=[0;U;2];
u=u(x)';
max(abs(u-U))

%%
clear all
clc
close all

f=@(x) (pi/2-x).^5.*cos(x)-10.*x;
a=0;
b=pi/2;
N=9;
h=(b-a)/N;
x=(a+b)/2;
fa=f(a);
fx=f(x);
for k=1:9
    if fa*fx<0
        b=x;
    else
        a=x;
    end
    x=(a+b)/2;
    fa=f(a);
    fx=f(x);
end
(a+b)/2

%%
clear all
clc
close all

x0=0;
xN=5;
N=100;
x=[0,1,2,3,4,5];
y=[1.2,2.3,1.9,2.1,0.8,3.3];
x0=linspace(x0,xN,N+1);
s=spline(x,y,x0);
h=(xN-x0)/N;
I=h/2*(s(1)+2*sum(s(2:N))+s(N+1))

%%
clear all
clc
close all

a=[2,1];
b=[4,2];
f=@(x,y) 3*atan((x-a(1)).^2+(y-a(2)).^2)+2*atan((x-b(1)).^2+(y-b(2)).^2);
x=linspace(0,5);
y=x;
[X,Y]=meshgrid(x,y);
surf(X,Y,f(X,Y))
hold on
plot3(1.86,0.92,f(1.86,0.92),'or')

%%
clear all
clc
close all

x0=0;
xN=4;
h=0.25;
N=(xN-x0)/h;
x0=linspace(x0,xN,N+1);
f=@(x,y) [ (2-2.*y(2)).*y(1); (y(1)-1).*y(2) ];
y0 = [0.25,0.25];
[t23,y23]=ode23(f, x0, y0);
[t45,y45]=ode45(f, x0, y0);
format long e
max(y45(:,1))
max(y23(:,1))
max(y45(:,2))
max(y23(:,2))
max(abs(y23(:,1)-y45(:,1)))
max(abs(y23(:,2)-y45(:,2)))


%%
clear all
clc
close all

A=diag(3*ones(1,8))+diag(-1*ones(1,6),2)+diag(-1*ones(1,4),4)+diag(-1*ones(1,2),6);
A=A+A';
b=ones(8,1);
x0=0;

D=tril(A);
C=A-D;
for k=1:200
    x=D\(b-C*x0);
    if (norm(x-x0)/norm(x)<1.0e-12)
        break;
    end
    x0=x;
end
k

%%
clear all
clc
close all

F = @(x) [x(1).^2+x(2).^2-1;
x(2).^2+x(3).^2-1;
x(3).^2+x(4).^2-1;
x(1)+x(2)-1;
x(2)+x(3)-1;
x(3)+x(4)-1];
for a = [5, 10]
fprintf('ITERAZIONE per a = %d', a);
x0 = a * ones(1, 4);
lsqnonlin(F, x0)
end

%%
clear all
clc
close all

f = @(t, y) [y(2); y(3); -5*y(3)-8*y(2)-4*y(1)+t^2+t+1];
t0 = 1;
tN = 6;
N = 38;
h = (tN -t0)/N;
t = linspace(t0, tN, N+1);
y0 = [3, 4, 5];
y = zeros(3, N+1);
y(:, 1) = y0;
for n = 1:N
y(:,n+1) = y(:,n) + h*f(t(n), y(:,n));
end
y

%%
clear all
clc
close all

f = @(x) exp(-x.^2);
a = -1;
b = 1;
I = integral(f, a, b);
for N = [4, 8, 16, 32]
x = linspace(a, b, N+1);
y = f(x);
h = (b -a)/N;
I_N = h/2 * (y(1) + 2*(sum(y(2:N))) + y(end));
N = 2*N;
x = linspace(a, b, N+1);
y = f(x);
h = (b -a)/(N);
I_2N = h/2 * (y(1) + 2*(sum(y(2:N))) + y(end));
abs(I_2N -I)/abs(I_N -I)
end

%%
clear all
clc
close all

F = @(x) [x(1).^2+x(2).^2-1;
x(2).^2+x(3).^2-1;
x(3).^2+x(4).^2-1;
x(1)+x(2);
x(2)+x(3);
x(3)+x(4)];
x0=4*ones(1,4);
lsqnonlin(F, x0)

%%
clc
close all
clear all

f=@(t,x)[x(2); 3*x(2)-10*x(1)+t];
y0=[pi;1];
t0=0;
tN=9;
h=0.05;
N=9/0.05;
x=linspace(t0,tN,N+1);
y=zeros(2,N+1);
y(:,1)=y0;
for n=1:N
    K1=f(x(n),y(:,n));
    K2=f(x(n+1),y(:,n)+h*K1);
    y(:,n+1)=y(:,n)+h/2*(K1+K2)
end
max(y(1,end))
    
%%
clear all
clc
close all

f=@(x,y) y-2*x+3;
x0=0;
y0=1;
N=5;
h=0.1;
xN=h*N+x0;
x=linspace(x0,xN,N+1);
y=zeros(1,N+1);
y(1)=y0;
for n=1:N
    K1=f(x(n),y(n));
    K2=f(x(n+1),y(n)+h*K1);
    y(n+1)=y(n)+h/2*(K1+K2);
end
y(end)
    
%%
clear all
clc
close all

f=@(x,y) -y+x.^2;
y0=sqrt(2);
x=20;
[x,y]=ode45(f,[0 20],sqrt(2));
x,y

%%
clear all
clc
close all

f=@(x,y)[y(2);3*y(2)-2*y(1)];
y0=[0;2];
N=100;
x0=0;
xN=1;
x=linspace(x0,xN,N+1);
h=1/N;
y=zeros(2,N+1);
y(:,1)=y0;
for n=1:N
    y(:,n+1)=y(:,n)+h*f(x(n),y(:,n));
end
y(1,end)
abs(2*exp(2)-2*exp(1)-y(1,end))

%%
clear all
clc
close all

f=@(x,y)[y(2);2.*y(2)-y(1)-x.^2];
y0=[pi;5];
x0=0;
xN=10;
h=0.1;
N=(xN-x0)/h;
x=linspace(x0,xN,N+1);
y=zeros(2,N+1);
y(:,1)=y0;
for n=1:N
    K1=f(x(n),y(:,n));
    K2=f(x(n+1),y(:,n)+h*K1);
    y(:,n+1)=y(:,n)+h/2*(K1+K2);
end
y(:,end)

%%
clear all
clc
close all

x0=-1;
xN=1;
N=1000;
h=(xN-x0)/N;
A=diag(6/h^2*ones(1,N-1))+diag(-3/h^2*ones(1,N-2),-1)+diag(-3/h^2*ones(1,N-2),1);
f=@(x)exp(pi*x);
x=linspace(x0,xN,N+1);
F=f(x(2:end-1))';
U=A\F;
max(U)

%%
clear all
clc
close all

f=@(x) sin(x);
x0=0;
xN=3;
N=100;
h=(xN-x0)/N;
A=diag(8/h^2*ones(1,N-1))+diag(-4/h^2*ones(1,N-2),-1)+diag(-4/h^2*ones(1,N-2),1);
x=linspace(x0,xN,N+1);
F=f(x(2:end-1))';
F(1)=F(1)+2/h^2;
F(N-1)=F(N-1)+4/(3*h^2);
U=A\F;
max(U)

%% DIFF FINITE - DERIVATA IN UN PUNTO CON DIFF FINITE
% 
clear all
close all
clc

x0=0;
xN=pi;
N=1000;
h=(xN-x0)/N;
A=diag(-6/h^2*ones(1,N-2),-1)+diag(12/h^2*ones(1,N-1))+diag(-6/h^2*ones(1,N-2),1);
A(1,1)=6/h^2;
x=linspace(x0,xN,N+1);
f=@(x)cos(x-log(x+2));
F=f(x(2:end-1))';
F(1)=F(1)-1/h;
F(end)=F(end)-6/h^2;
U=A\F;
max(U)

%%
clear all
close all
clc

x0=0;
xN=2*pi;
N=10;
h=(xN-x0)/N;
f=@(x) -sin(x);
A=diag((-2/h^2)*ones(1,N))+diag((1/h^2)*ones(1,N-1),-1)+diag((1/h^2)*ones(1,N-1),1);
A(end,end-1)=2/h^2;
A
x=linspace(x0,xN,N+1);
F=f(x(2:N+1))';
F(end)=F(end)-2/h;
U=A\F;

%%
clear all
clc
close all

N=100;
x0=0;
xN=1;
h=(xN-x0)/N;
A=diag((1/(5*h^2)+1/h)*ones(1,N-1))+diag((-1/(10*h^2)-1/h)*ones(1,N-2),-1)+diag((-1/(10*h^2))*ones(1,N-2),1);
x=linspace(x0,xN,N+1);
F=ones(1,N-1)';
U=A\F;
max(U)

%%
clear all
clc
close all

N=100;
x0=0;
xN=1;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
A=diag((10/h^2)*ones(1,N-1))+diag((-5/h^2)*ones(1,N-2),-1)+diag((-5/h^2)*ones(1,N-2),1);
f=@(x)log(3*x.^2+2);
F=f(x(2:end-1))';
U=A\F;
max(U)

%%
clear all
clc
close all

N=10;
x0=0;
xN=3;
h=(xN-x0)/N;
d=1/h^2*ones(N-1,1);
A=spdiags([-2*d 4*d -2*d], -1:1, N-1, N-1);
f=@(x)x.^(2/3);
x=linspace(x0,xN,N+1);
F=f(x(2:end-1))';
F(1)=F(1)+2/h^2;
F(end)=F(end)+4/h^2;
U=A\F;
max(U)

%%
clear all
clc
close all

x0=0;
xN=pi;
N=1000;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
d=1/h^2*ones(N,1);
A=spdiags([-2*d 4*d -2*d], -1:1, N, N);
A(end,end-1)=-4/h^2;
f=@(x) sin(x+cos(x));
F=f(x(2:end))';
F(1)=F(1)+2/h^2;
F(end)=F(end)-2/h;
U=A\F;
max(U)


%%
clear all
clc
close all

x0=0;
xN=2*pi;
N=10;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
d=1/h^2*ones(N,1);
A=spdiags([d -2*d d], -1:1, N, N);
A(end,end-1)=2/h^2;
f=@(x) cos(x);
F=f(x(2:end))';
F(1)=F(1)+1/h^2;
U=A\F


%%
clear all
clc
close all

x0=0;
xN=1;
N=100;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
d=ones(N-1,1);
A=spdiags([(-1/(1000*h^2)-1/h)*d (2/(1000*h^2)+1/h)*d (-1/(1000*h^2))*d], -1:1, N-1, N-1);
F=ones(N-1,1);
U=A\F;
max(U)

%%
clear all
clc
close all

N=1000;
x0=0;
xN=2;
h=(xN-x0)/N;
f=@(x) abs(cos(pi*x));
x=linspace(x0,xN,N+1);
d=1/h^2*ones(N-1,1);
A=spdiags([-2*d 4*d -2*d], -1:1, N-1, N-1);
F=f(x(2:end-1))';
U=A\F;
max(U)

%%
clear all
clc
close all

x0=0;
xN=pi;
N=1000;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
d=1/h^2*ones(N-1,1);
A=spdiags([-6*d 12*d -6*d], -1:1, N-1, N-1);
f=@(x) cos(sin(x));
F=f(x(2:end-1))';
F(1)=F(1)+3/h^2;
F(end)=F(end)+2/h^2;
U=A\F;
max(U)

%%
clear all
clc
close all

x0=0;
xN=pi;
N=1000;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
d=1/h^2*ones(N,1);
A=spdiags([-3*d 6*d -3*d], -1:1, N, N);
A(end,end-1)=-6/h^2;
f=@(x) x+cos(x);
F=f(x(2:end))';
F(1)=F(1)-3/h^2;
F(end)=F(end)+4/h;
U=A\F;
max(U)

%%
clear all
clc
close all

x0=0;
xN=2*pi;
N=10;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
d=1/h^2*ones(N,1);
A=spdiags([d -2*d d], -1:1, N, N);
A(end,end-1)=2/h^2;
f=@(x) cos(x);
F=f(x(2:end))';
F(1)=F(1)+1/h^2;
U=A\F;
U(5)

%%
clear all
clc
close all

x0=0;
xN=1;
N=100;
h=(xN-x0)/N;
d=ones(N-1,1);
A=spdiags([ (-1/(1000*h^2)-1/h)*d (2/(1000*h^2)+1/h)*d (-1/(1000*h^2))*d ], -1:1, N-1, N-1);
F=ones(N-1,1);
U=A\F;
max(U)

%%
clear all
clc
close all

n=256;
a=2;
d=ones(n,1);
A=spdiags([-d a*d -d], -1:1, n, n);
x=ones(1,n)';
b=A*x;
pcg(A,b,10^(-2),100)

%%
clear all
clc
close all

A=[6,-9,-2; -3, 6, 1; -6, -6, -6];
b=[3, 0, -7]';
x0=[0, 0, 0]';
D=diag(diag(A));
C=A-D;
for k=1:10
    x=D\(b-C*x0);
    x0=x;
end
x0

%%
clear all
clc
close all

A=[-9,10,-7;4,-2,-2;-2,3,-7];
b=[5,8,-3];
D=diag(diag(A));
C=A-D;
B=-inv(D)*C;
max(abs(eigs(B)))

%%
clear all
clc
close all

n=256;
a=2;
d=ones(n,1);
A=spdiags([d a*d d], -1:1, n, n);
A(1,n)=0.01;
A(n,1)=0.01;
x=ones(n,1);
b=A*x;
x0=[0,0,0]';
L=ichol(A);
pcg(A,b,[],100,L,L')

%%
clear all
clc
close all

n=150;
d=ones(n,1);
A=spdiags([-d 4*d -d], -1:1, n, n);
b=nthroot(1:1:n,3)';
x0=zeros(150,1);
D=diag(diag(A));
C=A-D;
for k=1:136
    x=D\(b-C*x0);
    if norm(x-x0,1)/norm(x,1)<10^(-7)
        break;
    end
    x0=x;
end
k

%%
clear all
close all
clc

n=256;
a=2;
d=ones(n,1);
A=spdiags([d a*d d], -1:1, n, n);
A(1,n)=0.01;
A(n,1)=0.01;
x=ones(n,1);
b=A*x;
L=ichol(A);
pcg(A,b)
pcg(A,b,[],[],L,L')

%%
clear all
clc
close all

A=[-1,9,6;3,-7,0;-8,-5,6];
b=[-2;-5;-10];
x0=zeros(3,1);
D=diag(diag(A));
C=A-D;
for k=1:9
    x=D\(b-C*x0);
    x0=x;
end
x0

%%
clear all
clc
close all

n=256;
a=2;
d=ones(n,1);
A=spdiags([-d a*d -d], -1:1, n, n);
b=ones(n,1);
pcg(A,b)

%%
clear all
clc
close all

A=[5,1,0;1,5,-3;0,1,8];
b=[6;3;9];
D=diag(diag(A));
C=A-D;
B=-inv(D)*C;
max(abs(eigs(B)))

%%
clear all
clc
close all

f=@(x) x.^3+sin(x);
a=-2;
b=4;
x0=2;
x1=-3;
for k=1:4
    x=x1-f(x1)*(x1-x0)/(f(x1)-f(x0));
    x0=x1;
    x1=x;
end
x

%%
clear all
close all
clc

f=@(x) x.^3+3*x+sin(pi*x);
df=@(x) 3*x.^2+3+pi*cos(pi*x);
a=-2;
b=4;
x0=2;
for k=1:6
    x=x0-f(x0)/df(x0);
    x0=x;
end
x0

%%
clear all
clc
close all

f=@(x)x.^3.*sin(x);
df=@(x) 3*x.^2.*sin(x)+x.^3.*cos(x);
x0=-3;
for k=1:5
    x=x0-f(x0)/df(x0);
    x0=x;
end
x0

%%
clear all
clc
close all

f=@(x) exp(x)-10;
x0=1;
df=@(x) exp(x);
for k=1:7
    x=x0-f(x0)/df(x0);
    x0=x;
end
x0

%%
clear all
clc
close all

f=@(x) 1./(x.^2+sin(x)+1);
x0=0;
kmax=10;
for k=1:kmax
    x=f(x0);
    e(k+1)=abs(x-x0);
    x0=x;
end
e(10)/e(9)
p=log(e(2:end))./log(e(1:end-1));
p'

%%
clear all
clc
close all

a=-1;
b=3;
f=@(x) x.^2+3*x-1;
x=(b+a)/2;
fx=f(x);
fa=f(a);
for k=1:6
    if fa*fx<0
        b=x;
    else
        a=x;
    end
    fa=f(a);
    x=(b+a)/2;
    fx=f(x);
end
x

%%
clear all
clc
close all

f=@(x) x.^3+3*x+sin(pi*x)+pi;
a=-2;
b=3;
fa=f(a);
x=(b+a)/2;
fx=f(x);
for k=1:11
    if fa*fx<0
        b=x;
    else
        a=x;
    end
    fa=f(a);
    x=(b+a)/2;
    fx=f(x);
end
x

%%
clear all
clc
close all

f=@(x) exp(-x.^2./3);
x0=1;
for k=1:5
    x=f(x0);
    x0=x;
end
x0

%%
clear all
clc
close all

f=@(x) (x(1)-2).^2+x(2).^2+sin(pi*x(1));
x0=[0;0];
fminunc(f,x0)

f=@(x,y) (x-2).^2+y.^2+sin(pi*x);
x0=linspace(-2,2);
y0=x0;
[X,Y]=meshgrid(x0,y0);
surf(X,Y,f(X,Y))

%%
clear all
clc
close all

f=@(x) x(1).^2+4*x(2).^2-x(1)*x(2)+5*x(3).^2+x(1)+x(2)-x(3);
df=@(x) [2*x(1)-x(2)+1; 8*x(2)-x(1)+1; 10*x(3)-1];
ddf=@(x) [2, -1, 0; -1, 8, 0; 0, 0, 10];
x0=[0;1;0];
for k=1:100
    x=x0-ddf(x0)\df(x0);
    if norm(x-x0)/norm(x)<=1.0e-07
        break;
    end
    x0=x;
end
k
f(x)

%%
clear all
clc
close all

f=@(x) abs(x(1))+abs(x(2))+abs(x(3))+abs(x(4))+abs(x(5));
x0=[1,1,1,1,1]';
A=[0,1,0,-1,0];
b=1;

[X,FVAL,EXITFLAG,OUTPUT]=fmincon(f,x0,A,b,[],[],[],[],@costraint);
X,FVAL,OUTPUT

% function [c,ceq] = costraint(x)
%     c=2-abs(x(1))+abs(x(3));
%     ceq=[];
% end

%%
clear all
clc
close all

f=@(x)exp(-1./(1+x.^2));
x0=-0.5;
[x,fval,exit,output]=fminsearch(f,x0)

%%
clear all
clc
close all

f=@(x) (x(1)-2).^2+x(2).^2+sin(pi*x(1));
x0=[-1,0];
[x,fval,exit,output]=fmincon(f,x0,[],[],[],[],[],[],@costraint)
f=@(x,y) (x-2).^2+y.^2+sin(pi*x);
[X,Y]=meshgrid(linspace(-1,1),linspace(-1,1));
surf(X,Y,f(X,Y))
hold on
plot3(x(1),x(2),f(x(1), x(2)),'or')

% function [c,ceq] = costraint(x)
%     c=x(1)+x(2).^2;
%     ceq=[];
% end

%%
clear all
clc
close all

f=@(x)10*sin(x(1)+x(2))+(x(1)-3).^2+(x(2)-2).^2-5;
x0=[0,0];
[x,f,exit,output]=fminsearch(f,x0)

%%
clear all
clc
close all

f=@(x) abs(x(1))+abs(x(2))+abs(x(3))+abs(x(4))+abs(x(5));
x0=[1,1,1,1,1]';
[x,f,exit,output]=fmincon(f,x0,[],[],[],[],[-Inf,-Inf,2,-Inf,-Inf],[Inf,Inf,Inf,Inf,Inf])
    
%%
clear all
clc
close all

f=@(x) x(1).^2+4*x(2).^2-x(1).*x(2)+5*x(3).^2+x(1)+x(2)-x(3);
x0=[0;1;0];
[x,f,exit,output]=fmincon(f,x0)

%%
clear all
close all
clc

f=@(x)exp(-1./(1+x.^2));
df = @(x) (2*exp(-1/(1+x.^2)).*x)./((1+x.^2).^2);
dff = @(x) (exp(-1./(1+x.^2)).*(2-6.*x.^4))./((1+x.^2).^4);
x0=-0.5;
for k=1:6
    x=x0-df(x0)./dff(x0);
    if abs(x-x0)<=1.0e-07
        break;
    end
    x0=x;
end
f(x)

%%
clear all
clc
close all

f=@(x) x(1).^2+(x(2)-1).^2+sin(pi*x(1))+sin(pi*x(2));
x0=[4,4]';
[x,val,exit,output]=fmincon(f,x0,[],[],[],[],[0,-Inf],[Inf,0],@costraint)

% function [c,ceq] = costraint(x)
%     c=x(1).*x(2)+1;
%     ceq=[];
% end

%%
clear all
clc
close all

f=@(x)10*sin(x(1)+x(2))+(x(1)-3).^2+(x(2)-2).^2-5;
x0=[5;2];
df=@(x)[10*cos(x(1)+x(2))+2*(x(1)-3); 10*cos(x(1)+x(2))+2*(x(2)-2)];
dff=@(x)[-10*sin(x(1)+x(2))+2, -10*sin(x(1)+x(2)); -10*sin(x(1)+x(2)), -10*sin(x(1)+x(2))+2];
for k=1:6
    x=x0-dff(x0)\df(x0);
    if norm(x-x0)/norm(x)<=1.0e-07
        break;
    end
    x0=x;
end
x0
k
f(x)
f=@(x1,x2)10*sin(x1+x2)+(x1-3).^2+(x2-2).^2-5;
[X,Y]=meshgrid(linspace(0,10),linspace(0,10));
surf(X,Y,f(X,Y))
hold on
plot3(x(1),x(2),f(x(1),x(2)),'or');
hold on
plot3(4.5887,3.5887,f(4.5887,3.5887),'or');

%%
clear all
clc
close all

f=@(x)x(1).^2+(x(2)-1).^2+sin(pi*x(1))+sin(pi*x(2));
x0=[0,0];
[x,f,exit,output]=fminunc(f,x0)

%%
clear all
clc
close all

f=@(x,y)sin(pi*x).*sin(2*pi*y);
x=linspace(-1,1,6);
y=x;
[X,Y]=meshgrid(x,y);
z1=interp2(X,Y,f(X,Y),0.7,0.3,'spline');
z2=interp2(X,Y,f(X,Y),0.7,0.3,'linear');
abs(z1-z2)

%%
clear all
clc
close all

f=@(x,y)exp(-(sin(2*pi*x)).^2-(sin(2*pi*y)).^2);
x=linspace(-1,1,9);
y=x;
[X,Y]=meshgrid(x,y);
x2=linspace(-1,1);
y2=x2;
[X2,Y2]=meshgrid(x2,y2);
z=interp2(X,Y,f(X,Y),X2,Y2,'linear');
figure(1)
surf(X2,Y2,z)
figure(2)
surf(X2,Y2,f(X2,Y2))

%%
clear all
clc
close all

f=@(x,y)exp(-(sin(pi*x)).^2-(sin(pi*y)).^2);
x=linspace(-1,1,5);
y=x;
[X,Y]=meshgrid(x,y);
x2=linspace(-1,1);
y2=x2;
[X2,Y2]=meshgrid(x2,y2);
z=interp2(X,Y,f(X,Y),X2,Y2,'linear');
figure(1)
surf(X2,Y2,z);
figure(2)
surf(X2,Y2,f(X2,Y2));

%%
clear all
clc
close all

f=@(x,y)sin(pi*x).*sin(2*pi*y);
x=linspace(-1,1,6);
y=x;
[X,Y]=meshgrid(x,y);
x2=linspace(-1,1);
y2=x2;
[X2,Y2]=meshgrid(x2,y2);
z=interp2(X,Y,f(X,Y),X2,Y2,'spline');
figure(1)
surf(X2,Y2,z);
figure(2)
surf(X2,Y2,f(X2,Y2));

%%
clear all
clc
close all
format long e
x=0:0.2:1;
y=@(x)(sqrt(x)-(x+1).^2)./(x.^2+3);
spline(x,y(x),0.48)

%%
clear all
clc
close all

f=@(x)cos(x);
x=linspace(0,pi,5);
polyval(polyfit(x,f(x),4),pi/7)

%%
clear all
clc
close all

x=[-5,0,2];
y=[9,2,2];
spline(x,[6,y,7],sqrt(0.4))

%%
clear all
clc
close all

f=@(x)exp(-1./x);
n=5;
nodi=6;
x=-cos((pi*(2*(0:1:6)+1))/(2*(nodi)));
x2=(2-1)/2*x+(2+1)/2;
polyfit(x2,f(x2),n)'

%%
clear all
clc
close all

N=28;
x0=2;
xN=3;
h=(xN-x0)/(2*N);
f=@(x)x.^4.*(3*x-3)-4*x.^2;
x=linspace(x0,xN,2*N+1);
I=h/3*(f(x(1))+4*sum(f(x(2:2:2*N)))+2*sum(f(x(3:2:2*N-1)))+f(x(2*N+1)))

%%
clear all
clc
close all

N=72;
x0=8;
xN=9;
h=(xN-x0)/N;
f=@(x)log(x)+x;
x=linspace(x0,xN,N+1);
I=h/2*(f(x(1))+2*sum(f(x(2:N)))+f(x(N+1)))

%%
clear all
clc
close all

N=6;
x0=-3;
xN=5;
h=(xN-x0)/(2*N);
x=linspace(x0,xN,2*N+1);
f=@(x)log(x.^4+x.^2+1);
y=f(x);
z1=h/3*(y(1)+4*sum(y(2:2:2*N))+2*sum(y(3:2:2*N-1))+y(2*N+1))
N=12;
x0=-3;
xN=5;
h=(xN-x0)/(2*N);
x=linspace(x0,xN,2*N+1);
f=@(x)log(x.^4+x.^2+1);
y=f(x);
z2=h/3*(y(1)+4*sum(y(2:2:2*N))+2*sum(y(3:2:2*N-1))+y(2*N+1))
abs(z1-z2)

%%
clear all
clc
close all
a=1;
for i=1:74
    a=a*(i^((-1)^(i+1)));
end
a

%%
clear all
clc
close all

f=@(x)exp(pi*x);
x0=-1;
xN=1;
N=1000;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
F=f(x(2:end-1))';
d=1/h^2*ones(N-1,1);
A=spdiags([-3*d 6*d -3*d],-1:1,N-1,N-1);
U=A\F;
max(U)

%%
clear all
clc
close all

x0=0;
xN=3;
N=100;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
d=1/h^2*ones(N-1,1);
A=spdiags([-2*d 4*d -2*d],-1:1,N-1,N-1);
f=@(x)x.^(2/3);
F=f(x(2:end-1))';
F(1)=F(1)+2/h^2;
F(end)=F(end)+4/h^2;
U=A\F;
max(U)

%%
clear all
clc
close all

x0=0;
xN=pi;
N=1000;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
f=@(x)x-3*cos(x);
F=f(x(1:end-1))';
F(end)=F(end)+1/(2*h^2);
d=-1/(2*h^2)*ones(N,1);
A=spdiags([d -2*d d],-1:1,N,N);
A(1,2)=-1/h^2;
U=A\F;
max(U)

%%
clear all
clc
close all

x0=0;
xN=2*pi;
N=10;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
d=1/h^2*ones(N,1);
A=spdiags([d -2*d d],-1:1,N,N);
A(end,end-1)=2/h^2;
f=@(x)-cos(x);
F=f(x(2:end))';
F(1)=F(1)-1/h^2;
U=A\F

%%
clear all
clc
close all

x0=0;
xN=1;
N=80;
h=(xN-x0)/N;
x=linspace(x0,xN,N+1);
d=ones(N-1,1);
A=spdiags([d*(-1/(10*h^2)-1/h),d*(1/(5*h^2)+1/h),d*(-1/(10*h^2))],-1:1,N-1,N-1);
F=-1*ones(N-1,1);
U=A\F;
min(U)

%%
clear all
clc
close all

f=@(x)abs(x(1)-1)+abs(x(2)+2)-1;
x0=[1;1];
A=[-1,1];
b=-4;
[x,f,exit,output]=fmincon(f,x0,A,b)
f=@(x1,x2)abs(x1-1)+abs(x2+2)-1;
x=linspace(-5,10);
y=x;
[X,Y]=meshgrid(x,y);
surf(X,Y,f(X,Y))
hold on
plot3(1.0714,-2.9286,2.0000e-08,'or')







    
    