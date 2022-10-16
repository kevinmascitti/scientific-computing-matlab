clear
close all
clc

n=linspace(-1,1,20);
f=@(x) exp(-x+x.^4+1);
c=polyfit(n,f(n),19);
z=linspace(-1,1,2000);
p=polyval(c,z);
max(p-f(z));

clear
close all
clc

M=magic(765);
I=10*eye(765);
A=M+I;
x=ones(765,1);
b=A*x;
x=A\b;
norm(b-A*x,inf);

clear
close all
clc

v=linspace(1,3,15);
w=linspace(0,5,15);
prod=0;
for i=1:15
prod=prod+v(i)*w(i);
end
prod

clear
close all
clc

H=hilb(7);
I=0.001*eye(7);
A=H+I;
x=[1 2 3 4 5 6 7]';
b=A*x;
xerr=A\b;
N=norm(x-xerr,inf);

clear
close all
clc

A=diag(6*ones(1,18))+diag(3*ones(1,17),1)+diag(3*ones(1,17),-1);
x=0;
for j=1:3
    b=linspace(0,j,18)';
    x=x+A\b;
end
norm(x,2);

clear
close all
clc

n=100;
for i=1:n
    for j=1:n
        a(i,j)=cos(1/min(i,j));
    end
end
cond(a,1)

clear
close all
clc

x=10.^(-4);
y=7-sqrt(49+x.^2);
y2=(-x.^2)./(7+sqrt(49+x.^2));
err=abs(y-y2)/abs(y2);

clear
close all
clc

n=linspace(-1,1,10);
f=@(x) exp(-x+x.^3);
c=polyfit(n,f(n),9);
z=linspace(-1,1,1500);
p=polyval(c,z);
max(f(z)-p);

clear
close all
clc

n=100;
for i=1:n
    for j=1:100
        a(i,j)=sin(min(i,j));
    end
end
cond(a,inf);

clear
close all
clc

A=diag(12*ones(1,20))+diag(4*ones(1,19),1)+diag(4*ones(1,19),-1);
x=0;
for j=1:4
    b=linspace(0,j,20)';
    x=x+A\b;
end
norm(x);

clear
close all
clc

A=magic(567)+10*eye(567);
x=ones(567,1);
b=A*x;
x=A\b;
norm(b-A*x,inf)

clear
close all
clc

x=10.^(-7);
y=2-sqrt(4+x.^2);
yesatto=(-x.^2)./(2+sqrt(4+x.^2));
err=abs(y-yesatto)/abs(yesatto);

clear
close all
clc

f=@(x) sin(3*cos(x)).*cos(3*sin(x));
z=linspace(-4,4,2000);
max(f(z));

clear
close all
clc

A=hilb(15)+0.01*eye(15);
x=[-1:-2:-29]';
b=A*x;
xgauss=A\b;
N=norm(x-xgauss,inf)



