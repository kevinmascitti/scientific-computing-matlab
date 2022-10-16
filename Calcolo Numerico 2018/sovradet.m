clear
close all
clc

t=[0,8,18];
v=[44,43,67];
A=[t' ones(3,1)];
b=v';
c_star=A\b;
p=polyval(c_star,2);

clear
close all
clc

x=[-2 -1.3 -1 -0.7 -0.4 -0.1];
y=[0.3 0.5 1.5 1.3 0.8 0.1];
A=[x' 2*ones(6,1)];
b=y';
cstar=A\b;
p=polyval(cstar,3);
abs(p-1.5)

clear
close all
clc

x=[0.34 0.19 0.25 0.61 0.47 0.35 0.83];
y=[0.58 0.54 0.91 0.28 0.75 1.17 0.38];
A=[x' 2*ones(7,1)];
b=y';
cstar=polyfit(x,y,2);

clear
close all
clc

xi=linspace(0,pi/2,30);
fx=@(x) x.*sin(x);
fxi=fx(xi);
c=polyfit(xi,fxi,2);

clear
close all
clc

x=[0,1,2];
y=[1,2,4];
A=[x' ones(3,1)];
b=y';
cstar=A\b;

clear
close all
clc

x=[0.2 1.14 0.54 0.87 1.25 2.36 0.19 0.54 0.51 0.33];
y=[1.25 2.36 2.58 1.87 2.68 3.41 0.65 0.47 1.36 1.25];
A=[x' ones(10,1)];
b=y';
c=A\b;

clear
close all
clc

M=pascal(15);
I=eye(15);
C=M+I;
A=C(:,1:10);
x=ones(10,1);
b=A*x;
[Q,R]=qr(A);
xqr=R\(Q'*b);
R=chol(A'*A);
yen=R'\(A'*b);
xen=R\yen;
eqr=norm(abs(xqr-x)/x,2);
een=norm(abs(xen-x)/x,2);

clear
close all
clc

M=magic(20);
I=eye(20);
B=M+I;
A=B(:,1:10);
x=ones(10,1);
b=A*x;
[Q,R]=qr(A);
y=Q\b;
xqr=R\y;
xbs=A\b;
eqr=norm(abs(xqr-x)/x,2);
ebs=norm(abs(xbs-x)/x,2);

clear
close all
clc

xi=linspace(0,5,22);
fx=@(x) (x.^2).*log(1+x);
fxi=fx(xi);
c=polyfit(xi,fxi,4);
p=polyval(c,2);

clear
close all
clc

A=[0.2 1.14 0.54 0.87 1.25 2.36 0.19 0.54 0.51 0.33; 1 1 1 1 1 1 1 1 1 1]';
b=[1.25 2.36 2.58 1.87 2.68 3.41 0.65 0.47 1.36 1.25]';
c=A\b;

clear
close all
clc

x=[-2 -1.3 -1 -0.7 -0.4 -0.1];
y=[0.3 0.5 1.5 1.3 0.8 0.1];
A=[x' ones(6,1)];
b=y';
c=polyfit(x,y,2);
abs(polyval(c,3)-1.5);

clear
close all
clc

M=pascal(15);
I=eye(15);
C=M+I;
A=C(:,1:10);
x=ones(10,1);
b=A*x;
[Q,R]=qr(A);
xqr=R\(Q'*b);
r=chol(A'*A);
yen=R'\(A'*b);
xen=R\yen;
eqr=norm(abs(xqr-x)/abs(x),2);
een=norm(abs(xen-x)/abs(x),2);

clear
close all
clc

xi=linspace(0,pi/2,30);
f=@(x) x.*sin(x);
yi=f(xi);
c=polyfit(xi,yi,2);

clear
close all
clc