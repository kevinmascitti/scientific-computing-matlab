% 9 es 8
clear all
close all
clc
format short e

f = @(x) x(1).^2-x(2).^2;
x0 = [2;2];
A = [-1 -1];
b = -1;

[x,fval,exit,output] = fmincon(f,x0,A,b,[],[],[],[],@vincolononlin)

function [c,ceq] = vincolononlin(x)
    c = [x(1).^2+x(2).^2-9;x(1).^2-x(2)];
    ceq = [];
end