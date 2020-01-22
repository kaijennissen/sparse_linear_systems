using LinearAlgebra


struct spMatrix
    p
    i
    x
end

# Algorithm 2.1 - sparse matrix time dense vector
function spM_x_deV(A, x)
    n = size(x, 1);
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
sparseM_times_denseVA,x) == [4.6 2.3 0; 0 1.5 0; 3.0 0 2.7] * x






# Algorithm 2.2 - sparse matrix times sparse matrix

function spM_x_spM(A, B)
    k = size(B.p,1) - 1; #maximum([maximum(A.i), maximum(B.i)])
    C = zeros(k, k);
    w = zeros(7);
    Cj = Vector{Int64}();
    x = zeros(k);
    #n = length(B.p)
    for j in 1:3
        #idxB = B.p[j]+1]:B.i[B.p[j+1]
        beta_lower = B.p[j];
        beta_upper = B.p[j+1];
        beta_range = collect((beta_lower+1):beta_upper)
        beta = B.i[beta_range];
        for t in beta
            b_idx = beta_range[findall(x -> x == t, beta)[]];
            bij = B.x[b_idx];
            if A.p[t+1] == A.p[t] # check that column has nonzeros
                continue
            else
                alpha_lower = A.p[t]+1;
                alpha_upper = A.p[t+1];
                alpha = collect((A.p[t]+1):A.p[t+1]);
            end
            for i in A.i[alpha_lower:alpha_upper]
                #a_index = findall(x -> x == i, A.i[alpha_lower:alpha_upper])[]
                b_index = beta[findall(x -> x == t, beta)[]];
                #row_ind = A.i[i]
                if w[i] < j
                    push!(Cj, i);
                    w[i] = j;
                end
                aij = A.x[A.p[t]+1]; #A.x[b_index]; #
                #bij = B.x[b_index];
                x[i] = x[i] + aij * bij;
            end
        end

        for i in Cj
            C[i, j] = x[i];
            x[i] = 0;
        end

    end

    C

end

A = spMatrix([0, 1, 2, 2, 3, 4, 5], [1, 2, 1, 1, 3], [1, 8, 4, 5, 18])
B = spMatrix([0, 2, 4, 5], [3, 4, 1, 5, 6], [7, 10, 2, 14, 18])

Juno.@enter spM_x_spM(A,B);



Adense = [1 0 0 4 5 0; 0 8 0 0 0 0; 0 0 0 0 0 18];
Bdense = [0 2 0; 0 0 0; 7 0 0; 10 0 0; 0 14 0; 0 0 18];
Adense*Bdense
