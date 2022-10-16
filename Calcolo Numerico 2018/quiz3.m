clear
close all
clc

A=hilb(10);
z=ones(1,10)';
w=z/norm(z);
p=2;
iter=4;
lambda=0;
[L,U,P]=lu(A-p*eye(10));
for i=1:iter
    y=L\(P*w);
    z=U\y;
    lambda=p+1/(w'*z);
    w=z/norm(z);
end
d=eigs(A,1,p);
err=abs(lambda-d)/abs(d)

clear
close all
clc

f=@(x) 1./((x.^4)+9);
z=linspace(1,2,5);
c=polyfit(z,f(z),4);

clear
close all
clc

x=linspace(0,1,10);
A=vander(x);
[U,S,V]=svd(A);
A7=U(:,1:7)*S(1:7,1:7)*V(:,1:7)';
norm=norm(A7,inf)

clear
close all
clc

H=hilb(10);
b=[2:2:20]';
x=H\b;

clear
close all
clc

f=@(x) exp(-x).*cos(x.^2);
z=linspace(-3,-2,1000);
plot(z,f(z));

clear
close all
clc

f=@(x) (x.^5).*log(x);
z=linspace(1,2,35);
fp=@(x) (5.*(x.^4)).*log(x)+x.^4;
x=linspace(1,2,500);
s=spline(z,[fp(z(1)) f(z) fp(z(35))],x);
max(abs(f(x)-s))