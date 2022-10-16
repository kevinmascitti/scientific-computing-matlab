clear
close all
clc

A=hilb(4);
[L,U,P]=lu(A)

clear
close all
clc

A=[4,6; 3/5,1];
cond(A,inf);

clear
close all
clc

A=diag(6*ones(18,1))+diag(3*ones(17,1),1)+diag(3*ones(17,1),-1);
y1=linspace(0,1,18);
y2=linspace(0,2,18);
y3=linspace(0,3,18);

x1=A\y1';
x2=A\y2';
x3=A\y3';

norm(x1+x2+x3,2);

clear
close all
clc

A=diag(4*ones(100,1))+diag(-1*ones(99,1),1)+diag(-1*ones(99,1),-1)+diag(-2*ones(90,1),10)+diag(-2*ones(90,1),-10);
cond(A,inf)
