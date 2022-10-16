clear
close all
clc

A=hilb(4);
[L,U,P]=lu(A)

clear
close all
clc

A=[4 6;3/5 1];
c=cond(A,inf)

clear
close all
clc

A=[6*pi 3 2 1; 3 7*pi 1 0; 2 1 6 0; 1 0 0 4];
c=chol(A)

clear
close all
clc

a=-1*ones(99,1);
b=-2*ones(99,1);
A=4*diag(ones(100,1),1)+diag(a,1)+diag(b,-1));
p=cond(A,inf)