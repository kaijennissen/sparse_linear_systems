using LinearAlgebra

# sparse Matrix struct
struct spMatrix
    p;
    i;
    x;

    function spMatrix(p, i, x)
        p = convert(AbstractArray{Int, 1}, p)
        i = convert(AbstractArray{Int, 1}, i)
        x = convert(AbstractArray{Float64, 1}, x)
        new(p, i, x)
    end
end


# convert from sparse to dense
function sparse2dense(X)
    n = length(X.p)-1;
    k = maximum(X.i);
    Y = zeros(k, n);
    for i in 1:n

        alpha_lower = X.p[i]+1;
        alpha_upper = X.p[i+1];

        if alpha_lower > alpha_upper
            continue
        end

        if alpha_upper > size(X.i, 1)
            alpha_upper = alpha_lower
        end

        #range=(X.p[i]+1):X.p[i+1]
        for j in collect(alpha_lower:alpha_upper);

            Y[X.i[j], i] = X.x[j]
        end
    end
    return Y
end


# convert from dense to sparse
function dense2sparse(X)
    n, k = size(X);
    X_nnz = X .!= 0;
    i = [];
    p = collect(Iterators.flatten(cumsum(sum(X_nnz, dims=1), dims=2)));
    pushfirst!(p, 0);
    x = X[X_nnz]

    for ii in 1:k
        append!(i, findall(x -> x == 1, X_nnz[:,ii]))
    end

    return spMatrix(p, i, x)
end


# Algorithm 2.1
"""
sparse matrix time dense vector
"""
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
@show spM_x_deV(A,x) == sparse2dense(A)*x;


# Algorithm 2.2
"""
sparse matrix times sparse matrix
"""
function spM_x_spM(A::spMatrix, B::spMatrix)
    n = (size(B.p, 1) - 1)::Int; #maximum([maximum(A.i), maximum(B.i)])
    #C = zeros(n, n)::AbstractArray{Float64, 2};
    w = zeros(6)::AbstractArray{Float64, 1};
    Cj = Vector{Int64}();
    Cp = zeros(n)::AbstractArray{Float64, 1};
    Cx = [];
    x = zeros(n)::AbstractArray{Float64, 1};

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
B = spMatrix([0, 2, 4, 5], [3, 4.0, 1, 5, 6], [7, 10.0, 2, 14, 18])
C = spM_x_spM(A,B);

@show sparse2dense(A)*sparse2dense(B) == sparse2dense(C);


# Algorithm 3.1
"""
lower triangular solve Lx=b with dense b
"""
function lu_solve(L, b)
    n = size(L.p, 1)-1;
    x = deepcopy(b);

    for j in 1:n
        alpha_lower = L.p[j]+2;
        alpha_upper = L.p[j+1];
        x_i = L.x[L.i[L.p[j]+1]];
        if alpha_lower > alpha_upper
            continue
        end

        if alpha_upper > size(L.i, 1)
            alpha_upper = alpha_lower
        end

        #alpha = L.i[]

        for i in alpha_lower:alpha_upper
            row = L.i[i];
            l_ij = L.x[i];
            x[row] = x[row] - l_ij*x[j];
        end
    end

    return x

end

L = [1 0 0; 2 1 0 ; -1 0 1];
x = [1; 2; 9];
b = L*x;
Lsp = dense2sparse(L);

@show lu_solve(Lsp, b) == x;


# Algorithm 3.2

"""
lower triangular solve Lx=b with sparse b
"""
function lu_solve_sparse(L, b)
    n = size(L.p, 1)-1;
    x = deepcopy(b);

    for j in 1:n
        alpha_lower = L.p[j]+2;
        alpha_upper = L.p[j+1];
        x_i = L.x[L.i[L.p[j]+1]];
        if alpha_lower > alpha_upper
            continue
        end

        if alpha_upper > size(L.i, 1)
            alpha_upper = alpha_lower
        end

        #alpha = L.i[]

        for i in alpha_lower:alpha_upper
            row = L.i[i];
            l_ij = L.x[i];
            x[row] = x[row] - l_ij*x[j];
        end
    end

    return x

end
