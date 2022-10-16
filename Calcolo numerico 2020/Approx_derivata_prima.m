%% approssimazione della derivata di f
clear all
close all
clc

f = @(x) exp(x);
df = @(x) exp(x);
x = 1;
k = 1:50;
h = 2.^(-k);
    % diminuendo la h, si avvicinano i
    % valori delle differenze finite con i 3 metodi
    % e gli errori relativi tendono a 0
    % h si dimezza e gli errori si dimezzano
    % diff finite servono a far convergere l'err_ass a 0 linearmente
% diff in avanti
Da = (f(x+h)-f(x))./h;
% diff all'indietro
Di = (f(x)-f(x-h))./h;
% diff centrata
Dc = (f(x+h)-f(x-h))./(2*h);

err_Da = abs(df(x)-Da);
err_Di = abs(df(x)-Di);
err_Dc = abs(df(x)-Dc);
disp('    h    Diff_avanti Diff_indietro Diff_centrale')
[h' err_Da' err_Di' err_Dc']

loglog(h,err_Da,'r',h,err_Di,'b',h,err_Dc,'m','linewidth',2)
legend('Diff-avanti', 'Diff-indietro', 'Diff-centrale')
% se k aumenta, l'errore diminuisce
% le formule di differenze finite sono soggette a
% cancellazione numerica perché c'è una differenza se aumento troppo k
