function lambda=potenze(A,z,toll,n_max)
Deltalambda=inf;
lambda=0;
n=1;
while (n<=n_max&Deltalambda>toll)
    w=z/norm(z);
    z=A*w;
    lambda(n+1)=w*z;
    Deltalambda=(lambda(n+1)/lambda(n))/abs(lambda(n+1));
    n=n+1;
end