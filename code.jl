using LinearAlgebra


struct spMatrix
    p
    i
    x
end


Adense = [4.6 2.3 0; 0 1.5 0; 3.0 0 2.7]
x = [1; 3; 4]
Adense*x

A = spMatrix([0, 2, 4, 5], [1, 3, 1, 2, 3], [4.6, 3.0, 2.3, 1.5, 2.7])

# Algorithm 2.1 - sparse matrix time dense vector
y = zeros(3);
n = size(x,1);
nnz = length(A.x);
for i in 1:n
    for j in (A.p[i]+1):(A.p[i+1])
        row = A.i[j]
        y[row] = y[row] + A.x[j]*x[i];
    end
end
