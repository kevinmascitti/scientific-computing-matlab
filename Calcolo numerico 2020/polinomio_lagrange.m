function li = polinomio_lagrange(x,i,z) 
% x rappresenta il vettore di valori
% i rappresenta l'indice fisso e ind varia
% z rappresenta il punto in cui valuto il polinomio
li=1;
for ind=1:length(x)
    if ne(ind,i)
    li=li.*(z-x(ind))./(x(i)-x(ind));
    end
end
end