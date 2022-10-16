function ex = esponenziale(x,toll)
ex=0;
i=0;
while (x^i/factorial(i)) >= toll
    ex = ex + x^i/factorial(i);
    i = i+1;
end

err_rel = abs(ex -exp(x))/abs(exp(x))
% sviluppo migliore per numero vicino allo 0