clear
close all
clc

A=rand(100);
b=sum(A,2);

tic
[L,U,P]=lu(A);
for i=1:30
    x=U\(L\(P*b));
    b=x;
end
TempoLU=toc

b=sum(A,2);

tic
for i=1:30
    x=A\b;
    b=x;
end
TempoMAT=toc