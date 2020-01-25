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

function spM_x_spM(A::spMatrix, B::spMatrix)
    n = size(B.p, 1) - 1; #maximum([maximum(A.i), maximum(B.i)])
    C = zeros(n, n);
    w = zeros(6);
    Cj = Vector{Int64}();
    Cp = zeros(n);
    Cx = [];
    x = zeros(n);

    for j in 1:n

        beta_lower = B.p[j];
        beta_upper = B.p[j+1];
        beta_range = collect((beta_lower+1):beta_upper)
        beta = B.i[beta_range];

        for t in beta
            b_idx = beta_range[findall(x -> x == t, beta)[]];
            bij = B.x[b_idx];

            alpha_lower = A.p[t]+1;
            alpha_upper = A.p[t+1];

            if alpha_lower > alpha_upper
                continue
            end

            if alpha_upper > size(A.i,1)
                alpha_upper = alpha_lower
            end

            alpha_range = collect(alpha_lower:alpha_upper)
            global alpha = A.i[alpha_range];

            for i in alpha

                a_index = alpha_range[findall(x -> x == i, alpha)[]];

                if w[i] < j
                    push!(Cj, i);
                    w[i] = j;
                    Cp[j] += 1
                end

                aij = A.x[a_index];
                x[i] = x[i] + aij * bij;

            end
        end

        push!(Cx, deepcopy(x[alpha]));
        x = zeros(n);

    end

    push!(w, 0)
    Cp = cumsum([0; Cp])
    Cx = collect(Iterators.flatten(Cx))
    return spMatrix(Cp, Cj, Cx)
end

A = spMatrix([0, 1, 2, 2, 4, 5, 7], [1, 2, 1, 3, 1, 3], [1, 8, 4, 16, 5, 18])
B = spMatrix([0, 2, 4, 5], [3, 4, 1, 5, 6], [7, 10, 2, 14, 18])
C = spM_x_spM(A,B);

Juno.@enter spM_x_spM(A,B);



Adense = [1 0 0 4 5 0; 0 8 0 0 0 0; 0 0 0 16 0 18];
Bdense = [0 2 0; 0 0 0; 7 0 0; 10 0 0; 0 14 0; 0 0 18];
Adense*Bdense
