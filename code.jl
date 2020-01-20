using LinearAlgebra


struct spMatrix
    p
    i
    x
end

# Algorithm 2.1 - sparse matrix time dense vector
function sparse_mult(A, x)
    n = size(x,1);
    y = zeros(n);
    nnz = length(A.x);
    for i in 1:n
        for j in (A.p[i]+1):(A.p[i+1])
            row = A.i[j]
            y[row] = y[row] + A.x[j]*x[i];
        end
    end
    return y
end

A = spMatrix([0, 2, 4, 5], [1, 3, 1, 2, 3], [4.6, 3.0, 2.3, 1.5, 2.7])
x = [1; 3; 5]
sparse_mult(A,x) == [4.6 2.3 0; 0 1.5 0; 3.0 0 2.7] * x






# Algorithm 2.2 - 
