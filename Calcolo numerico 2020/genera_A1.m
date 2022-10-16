function A = genera_A1(n)
A = zeros(n);
for i = 1:n;
    for j = 1:n
        A(i,j) = max(i,j);
    end
end

