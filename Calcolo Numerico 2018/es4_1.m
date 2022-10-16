clear
close all
clc

n=100;
p=n:-1:1;

A1_1=ones(n,1)*p;
A1_2=A1_1-diag(ones(n-1,1),-1);
A1=triu(A1_2,-1);
A2=triu(A1_2)+triu(A1,1);

[S1]=eig(A1);
A1pert=A1;
A1pert(end,:)=A1pert(end,:)+1e-10;
S1pert=eig(A1pert);

figure
hold on
grid on %visualizza griglia sullo schermo
plot(real(S1),imag(S1),'bo') %real prende i valori reali e imag quelli immaginari del vettore S1
plot(real(S1pert),imag(S1pert),'r*')

[S2]=eig(A2);
A2pert=A2;
A2pert(end,:)=A2pert(end,:)+1e-10;
S2pert=eig(A2pert);

figure
hold on
grid on
plot(real(S2),imag(S2),'bo')
plot(real(S2pert),imag(S2pert),'r*')