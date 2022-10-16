function [y] = serie_armonica(m)
val_esatto = pi^2/6;
t = 4;
y = 0;
for k=1:m
    term_gen = chop(1/chop(k^2, t), t); % arrotondo k^2 a t cifre significative della mantissa
    y = chop(y + term_gen, t);
end
y
err_rel = abs(y-val_esatto)/abs(val_esatto)

y2=0;
for k=m:-1:1
    term_gen = chop(1/chop(k^2, t), t); % arrotondo k^2 a t cifre significative della mantissa
    y2 = chop(y2 + term_gen, t);
end
y2
err_rel2 = abs(y2-val_esatto)/abs(val_esatto)