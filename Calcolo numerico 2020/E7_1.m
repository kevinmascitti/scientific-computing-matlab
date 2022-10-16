% 7 es 1
clear all
close all
clc

% a = -5;
% b = 5;
% f = @(x) 1./(1+x.^2);
% I_es = 2*atan(5);

a = -1;
b = 1;
f = @(x) sqrt(1-x.^2);
I_es = pi/2;

for k = 1:12   
    m = 2^k;
    vm(k) = m;
    I_trap(k) = trapezi(f,a,b,m);
    err_rel_trap(k) = abs(I_es - I_trap(k))/abs(I_es);
    I_simp(k) = simpson(f,a,b,m);
    err_rel_simp(k) = abs(I_es - I_simp(k))/abs(I_es);
end

% ordine sperimentale di convergenza
EOC_trap = log2(err_rel_trap(1:end-1)./err_rel_trap(2:end));
EOC_simp = log2(err_rel_simp(1:end-1)./err_rel_simp(2:end));

disp('val. funz.   err_rel      val.fun.   err_rel')
[(vm+1)' err_rel_trap' (2*vm+1)' err_rel_simp']

disp('EOC  trap     EOC simpson')
[EOC_trap' EOC_simp']