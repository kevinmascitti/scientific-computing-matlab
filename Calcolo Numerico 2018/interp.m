clear
close all
clc

x=linspace(1,2,25);
y=log(x);
z=linspace(1,2,1000);
s=spline(x,y,z);

err=norm(log(z)-s,inf)/norm(log(z),inf)

clear
close all
clc

x=linspace(0,2*pi);
z=linspace(0,2*pi,8);
y=3.*z.*cos(z).*sin(z);
c=polyfit(z,y,7);
val=polyval(c,pi/2)

clear
close all
clc

x=linspace(0,pi,5);
y=sin(x);
c=polyfit(x,y,4);
z=polyval(c,pi/8)

clear
close all
clc

x=linspace(-1,1,1000);
y=3*x.*exp(4*x.^2).*cos(2*pi*x/3)+50*(1+sin(7*x));

xi=linspace(-1,1,7)
yi=3*xi.*exp(4*xi.^2).*cos(2*pi*xi/3)+50*(1+sin(7*xi))
c=polyfit(xi,yi,6)
z=polyval(c,xi)

plot(x,y,'r',xi,z,'b')

clear
close all
clc

x=[0.0, 0.5, 1.0, 1.5, 2.0];
y=(sin(x)-(x+1).^2)./(x.^2+3);
s=spline(x,y,1.97)

clear
close all
clc

x=linspace(0,1,1000);
y=atan(x.*(x+1));
xi=linspace(0,1,8);
yi=atan(xi.*(xi+1));
c=polyfit(xi,yi,7);
z1=polyval(c,0.5);
z2=polyval(c,0.7);
err1=atan(0.5*(0.5+1))-z1
err2=atan(0.7*(0.7+1))-z2

clear
close all
clc

x=linspace(-1,1,100);
y=exp(x)./(x.^2+1);
xi=linspace(-1,1,7);
yi=exp(xi)./(xi.^2+1);
c=polyfit(xi,yi,6);
p=polyval(c,xi);
z=polyval(c,x);
err=norm(y-z,inf)

clear
close all
clc

x=[-5,4,5,11];
y=[6,2,4,10];
s=spline(x,[10 y 4], sqrt(1.8))

clear
close all
clc

xi=linspace(0,2*pi,4);
yi=cos(xi);
c=polyfit(xi,yi,3);
p=polyval(c,pi/8)

clear
close all
clc
%SPLINE CUBICA VINCOLATA
n=25;
f=@(x) x.^3.*cos(x);
x=linspace(0,1,n);
y=f(x);
d=@(x) 3.*x.^2.*cos(x)-x.^3.*sin(x);
z=linspace(0,1,1000);
s=spline(x,[d(0) y d(1)],z);
err=max(abs(f(z)-s));

clear
close all
clc