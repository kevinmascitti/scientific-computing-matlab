clear
close all
clc

n=100;

for i=1:n
    for j=1:n
        A(i,j)=i*max(i,j);
    end
end

[L,U,P]=lu(A);
invLU=inv(U)*inv(L)*P;
invMAT=inv(A);

errrel=norm(invLU-invMAT,inf)\norm(invMAT,inf)