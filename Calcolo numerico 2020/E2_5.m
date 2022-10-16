% 2 es 5
clear all
close all
clc
n = 50;
A = eye(n);
A(1,:) = 1;
A(:,1) = 1;
[L,U,P] = lu(A);
subplot(1,4,1)
spy(A)
title('spy(A)')
subplot(1,4,2)
spy(L)
title('spy(L)')
subplot(1,4,3)
spy(U)
title('spy(U)')
subplot(1,4,4)
spy(P)
title('spy(P)')