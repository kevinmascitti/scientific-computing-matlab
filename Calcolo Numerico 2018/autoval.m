
clear
close all
clc

for i=1:12
for j=1:12
if i==j
    a(i,j)=2*i;
end
if i>j
    a(i,j)=-2/j;
end
if i<j
    a(i,j)=2/j;
end
j=j+1;
end
i=i+1;
end
x=max(abs(eig(a)));

clear
close all
clc

A=hilb(18);
iter=6;

for m=1:iter
    [Q,R]=qr(A);
    A=R*Q;
end
d=diag(A);

lambda=eig(A);

max(abs(lambda-d));

clear
close all
clc

A=hilb(6);
iter=4;
z=ones(1,6)';

n=size(A);
w=z/norm(z);
lambda=0.2;
[L,U,P]=lu(A-0.2*eye(n));
for m= 1:iter
    y=L\(P*w);
    z=(U\y);
    lamp=0.2+1/(w'*z);
    w=z/norm(z);
lambda=lamp;
end
[X,D]=eigs(A,1,0.2);

err=abs(diag(D)-lambda)/abs(diag(D));

clear
close all
clc

A=pascal(8);
[U,S,V]=svd(A);
A5=U(:,1:5)*S(1:5,1:5)*V(:,1:5)';
S;

clear
close all
clc

A=diag(9*ones(1,24))+diag(2*ones(1,23),+1)+diag(-2*ones(1,23),-1);
b=linspace(0,1,24)';
[U,S,V]=svd(A);
y=S\(U'*b);
x=V*y;
norm(x,2)+norm(y,2)

clear
close all
clc

for i=1:12
    for j=1:12
        if i==j
            A(i,j)=2*i;
        end
        if i<j
            A(i,j)=-2/j;
        end
        if i>j
            A(i,j)=2/j;
        end
        j=j+1;
    end
    i=i+1;
end

x=eigs(A,1);

clear
close all
clc

for i=1:100
    for j=1:100
        A(i,j)=(1+i+j)/(1+abs(i-j));
    end
end

%METODO DELLE POTENZE INVERSE PER AUTOVALORE PIU VICINO A P
[X,D]=eigs(A,1,100);
norm(X,1);

clear
close all
clc

%METODO DELLE POTENZE INVERSE
n=6;
p=0.2;
A=hilb(n);
iter=4;
lambda=p;
z=ones(1,n)';
w=z/norm(z);
[L,U,P]=lu(A-lambda*eye(n));
for i=1:iter %Az=w PAz=Pw LUz=Pw
    y=L\(P*w);
    z=U\y;
    lambdap=p+1/(w'*z);
    w=z/norm(z);
    lambda=lambdap;
end

real=eigs(A,1,p);
err=abs(real-lambda)/abs(real);

clear
close all
clc

A=[2,3,4;3,4,6;1,0,3];
B=A'*A;

clear
close all
clc

%DISTANZA MATRICE DA MATRICI DI RANGO K
A=pascal(8);
[U,S,V]=svd(A);
A5=U(:,1:5)*S(1:5,1:5)*V(:,1:5)';
norm(A-A5,2);

clear
close all
clc

%SVD RISOLUZIONE SISTEMA
A=diag(5*ones(1,16))+diag(7*ones(1,15),1)+diag(-7*ones(1,15),-1);
b=linspace(0,1,16)';
[U,S,V]=svd(A);
y=S\(U'*b);
x=V*y;
norm(x)+norm(y);

clear
close all
clc

A=hilb(12);
iter=8;
for i=1:iter
    [Q,R]=qr(A);
    A=R*Q;
end

clear
close all
clc

A=hilb(12);
z=ones(1,12)';
w=z/norm(z);
lambda=0;
iter=7;
for i=1:iter
    z=A*w;
    lambdamax=w'*z;
    w=z/norm(z);
    lambda=lambdamax;
end
err=abs(eigs(A,1)-lambda)/abs(eigs(A,1));

clear
close all
clc

x=linspace(0,1,10);
A=vander(x);
[U,S,V]=svd(A);
A7=U(:,1:7)*S(1:7,1:7)*V(:,1:7)';
norm(A7,inf);

clear
close all
clc

A=diag(-4*ones(1,30))+diag(3*ones(1,29),1)+diag(-3*ones(1,29),-1);
b=linspace(0,1,30)';
[U,S,V]=svd(A);
y=S\(U'*b);
x=V*y;
norm(x)+norm(y);

clear
close all
clc

iter=4;
p=2;
A=pascal(6);
z=ones(1,6)';
w=z/norm(z);
lambda=p;
[L,U,P]=lu(A-lambda*eye(6));
for i=1:iter
    y=L\(P*w);
    x=U\y;
    lambdap=p+1/(w'*z);
    w=z/norm(z);
    lambda=lambdap;
end

real=eigs(A,1,p);
err=abs(lambda-real)/abs(real)

clear
close all
clc

A=hilb(18);

for i=1:6
    [Q,R]=qr(A);
    A=R*Q;
end

d=diag(A);
lambda=eig(A);
max(abs(lambda-d))

clear
close all
clc

n=50;
A=hilb(n);
condeig(A)

clear
close all
clc

x=linspace(0,1,10);
A=vander(x);
[U,S,V]=svd(A);
A7=U(:,1:7)*S(1:7,1:7)*V(:,1:7)';
norm(A7,inf)

clear
close all
clc

A=diag(9*ones(1,24))+diag(2*ones(1,23),1)+diag(-2*ones(1,23),-1);
b=linspace(0,1,24)';
[U,S,V]=svd(A);
y=S\(U'*b);
x=V*y;
norm(x)+norm(y)

clear
close all
clc

iter=4;
p=2;
A=pascal(6);
z=ones(1,6)';
w=z/norm(z);
lambda=p;
[L,U,P]=lu(A-lambda*eye(6));
for i=1:iter
    y=L\(P*w);
    z=U\y;
    lambdap=p+1/(w'*z);
    w=z/norm(z);
    lambda=lambdap;
end

real=eigs(A,1,p);
err=abs(lambda-real)/abs(real)

clear
close all
clc

for i=1:12
    for j=1:12
        if i==j
            a(i,j)=2*i;
        end
        if i<j
            a(i,j)=-2/j;
        end
        if i>j
            a(i,j)=2/j;
        end
    end
end

max(eig(a))

clear
close all
clc

iter=24;
x=linspace(-1,1,10);
A=vander(x);
z=ones(1,10)';
w=z/norm(z);
[X,D]=eigs(A,1)

clear
close all
clc
%METODO QR
iter=6;
A=hilb(18);
for i=1:iter
    [Q,R]=qr(A);
    A=R*Q;
end
lambda=eig(A);
d=diag(A);
err=max(abs(lambda-d))

clear
close all
clc

n=6;
iter=4;
p=0.2;
lambda=0;
A=hilb(n);
z=ones(1,n)';
w=z/norm(z);
[L,U,P]=lu(A-p*eye(n));
for i=1:iter
    y=L\(P*w);
    z=U\y;
    lambdap=p+1/(w'*z);
    lambda=lambdap;
    w=z/norm(z);
end
real=eigs(A,1,p);
err=abs(real-lambda)/abs(real)

clear
close all
clc

%METODO DELLE POTENZE
iter=24;
x=linspace(-1,1,10);
A=vander(x);
z=ones(1,10)';
w=z/norm(z);
lambda=0;
for i=1:iter
    z=A*w;
    lambda=w'*z;
    w=z/norm(z);
end
w(3)

n=55;
A=hilb(n);
condeig(A)

clear
close all
clc

A=pascal(8);
[U,S,V]=svd(A);
A5=U(:,1:5)*S(1:5,1:5)*V(:,1:5)';
norm(A-A5)

clear
close all
clc

A=diag(6*ones(1,18))+diag(3*ones(1,17),1)+diag(-3*ones(1,17),-1);
b=linspace(0,1,18)';
[U,S,V]=svd(A);
%A=USV' Ax=USV'x=b SV'x=U'b Sy=U'b
y=S\(U'*b);
x=V*y;
norm(x)+norm(y)

clear
close all
clc

n=12;
A=hilb(n);
z=ones(n,1);
w=z/norm(z);
lambda=0;
for i=1:7
z=A*w;
lambda=w'*z;
w=z/norm(z);
end
aval=eig(A);
m=max(abs(aval));
err=abs(m-lambda)/abs(m);

clear
close all
clc


