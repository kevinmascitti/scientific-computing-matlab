clear
close all
clc

tol=1e-10;
z=[1 2 3]';
itermax=100;

A1=[1 2 0; 1 0 0; 0 1 0];
lambda=potenze(A1,z,tol,itermax);
lambdaref=max(abs(eig(A1)))
(abs(lambdaref-lambda(end))/abs(lambdaref))