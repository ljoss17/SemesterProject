#etprecision(700)

log_sum_table = []

function fill_log_sum_table(C)
    old = 0.0
    for i in 1:C
        push!(log_sum_table, log(i)+old)
        old = log_sum_table[i]
    end
end

function compute_p(G, N)
    if G > N || N < 1
        throw(ArgumentError("Set G must be smaller than N"))
    end
    p = 1-(1 - G/N)^2
    return p
end

function sum_body(k, p, C)
    res = binomial(BigInt(C), BigInt(k))*(1-p)^(k*(C-k))
    return res
end

function log_binomial_coefficient(n, k)
    if n == k || k == 0
        return log(1)
    end
    if length(log_sum_table) < n+1
        fill_log_sum_table(n+1)
    end
    #println("test2 : ", length(log_sum_table))
    #println("n : ", n, " k : ", k)
    return log_sum_table[n] - log_sum_table[k] - log_sum_table[n-k]
end

function log_sum_body(k, p, C)
    res = log_binomial_coefficient(C, k) + (log(1-p) * (k*(C-k)))
    return Float64(exp(res))
end

function murmur_totality(G, C, N::Int64=524288, f::Float64=0.2)
    if f < 0 || f > 1
        throw(ArgumentError("Fraction of Byzantine processes f must be between 0 and 1"))
    end
    P = compute_p(G, N)
    ϵ = Float64(0.0)
    for k in 1:ceil(Int, C/2)
        ϵ = ϵ + log_sum_body(k, P, C)
    end
    return ϵ
end

function run_murmur(G, N::Int64=524290, f::Float64=0.2)
    C = floor(Int, (1-f)*N)
    fill_log_sum_table(C)
    for n in 1:G
        ϵ = murmur_totality(n, C)
        println("(", n, ", ", ϵ, ")")
    end
end
