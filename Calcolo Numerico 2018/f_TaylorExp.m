function u = f_TaylorExp(x,tolleranza)
u=0;
i=0;
termine=1;
while(i<=16 & termine>=tolleranza)
    u=u+termine;
    i=i+1;
    termine=x.^i/factorial(i);
end