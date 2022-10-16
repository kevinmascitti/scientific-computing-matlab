clear
close all
clc

f=@(x) cos(x);
n=linspace(-1,1,12);
z=linspace(-1,1,120);
s=spline(n,f(n),z);
max(abs(f(z)-s));

clear
close all
clc

x=10.^(-9);
f=sqrt((exp(x)-1)./x);
y=sqrt((x+(x.^2)./2)./x);
err=abs(f-y)/abs(y);

clear
close all
clc

f=@(x) x.*cos(x-2*pi);
z=linspace(3,6,9);
c=polyfit(z,f(z),2);
polyval(c,5.8)

clear
close all
clc

f=@(x) cos(1./x);
z=linspace(0.2,1,3);
abs(f(0.8)-polyval(polyfit(z,f(z),2),0.8));

clear
close all
clc

A=diag(12*ones(1,20))+diag(4*ones(1,19),1)+diag(4*ones(1,19),-1);
b=linspace(0,1,20)';
R=chol(A);
%R'*Rx=b
y=R'\b;
x=R\y;
norm(x+y,inf)

clear
close all
clc

x=linspace(-1,1);
f=exp(-x.^2).*sin(5.*x);
c=roots(f);

clear
close all
clc

A=hilb(100);
k=0;
for i=1:100
    for j=1:100
        if A(i,j)<0.03
            k=k+A(i,j);
        end
    end
end
k

clear
close all
clc

f=@(x) x.^4-2*x.^2+1;
z=linspace(0,1,8);
polyfit(z,f(z),7)

clear
close all
clc

f=@(x) exp(-cos(5.*x)).*x-1;
z=linspace(0,2,1000)
plot(z,f(z));