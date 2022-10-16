clear
close all
clc

n=5;
for i=1:n
    for j=1:n
        A(i,j)=max(i,j);
    end
end
b=sum(A,2);
[L,U,P]=lu(A);
y=L\(P*b);  %Ly=Pb
x=U\y;      %Ux=y

xref=ones(size(x));
errrel=norm(x-xref,inf)/norm(xref,inf)

