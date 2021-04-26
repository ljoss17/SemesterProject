log_sum_table = []

function fill_log_sum_table(C)
    # Fill table of logarithmic value of factorial i
    old = 0.0
    for i in 1:C
        push!(log_sum_table, log(i)+old)
        old = log_sum_table[i]
    end
end

function compute_p(G, N)
    # Probability p of any two correct processes being connected in a
    # Erdos–Rényi graph
    if G > N || N < 1
        throw(ArgumentError("Set G must be smaller than N"))
    end
    p = 1-(1 - G/N)^2
    return p
end

function log_binomial_coefficient(n, k)
    # Logarithmic value of a binomial coefficient of n and k.
    if n == k || k == 0
        return log(1)
    end
    if length(log_sum_table) < n+1
        fill_log_sum_table(n+1)
    end
    return log_sum_table[n] - log_sum_table[k] - log_sum_table[n-k]
end

function log_sum_body(k, p, C)
    res = log_binomial_coefficient(C, k) + (log(1-p) * (k*(C-k)))
    return Float64(exp(res))
end

function murmur_totality(G, N::Int64=524288, f::Float64=0.2)
    if f < 0 || f > 1
        throw(ArgumentError("Fraction of Byzantine processes f must be between 0 and 1"))
    end
    C = floor(Int, (1-f)*N)
    fill_log_sum_table(C)
    P = compute_p(G, N)
    ϵ = Float64(0.0)
    for k in 1:ceil(Int, C/2)
        ϵ = ϵ + log_sum_body(k, P, C)
    end
    return ϵ
end

function run_murmur(G, N::Int64=524290, f::Float64=0.2)
    for n in 1:G
        ϵ = murmur_totality(n)
        println("(", n, ", ", ϵ, ")")
    end
end
