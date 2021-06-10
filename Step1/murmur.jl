mutable struct SumTable
    table::Array{Float128,1}
    highest::Int64
end

log_sum_table = SumTable(zeros(1), 0)

function kahan_summation(values)::Float128
    # Kahan summation algorithm on array of values
    res::Float128 = 0.0
    c::Float128 = 0.0
    y::Float128 = 0.0
    t::Float128 = 0.0
    for i in values
        y = i - c
        t = res + y
        c = (t-res)-y
        res = t
    end
    return res
end

function fill_log_sum_table(C)
    # Fill table of logarithmic value of factorial i
    log_sum_table.table = zeros(C)
    old::Float128 = 0.0
    for i in 1:C
        log_sum_table.table[i] = log(i)+old
        old = log_sum_table.table[i]
    end
    log_sum_table.highest = C
end

function compute_p(G, N)::Float128
    # Probability p of any two correct processes being connected in a
    # Erdos–Rényi graph
    if G > N || N < 1
        throw(ArgumentError("Set G must be smaller than N"))
    end
    p::Float128 = 1-(1 - G/N)^2
    return p
end

function log_binomial_coefficient(n, k)::Float128
    # Logarithmic value of a binomial coefficient of n and k.
    if n == k || k == 0
        return log(1)
    end
    if log_sum_table.highest < n
        fill_log_sum_table(n)
    end
    return log_sum_table.table[n] - log_sum_table.table[k] - log_sum_table.table[n-k]
end

function log_sum_body(k, p, C)::Float128
    # Compute the body of the sum for Murmur totality using log and exp.
    res::Float128 = log_binomial_coefficient(C, k) + (log(1-p) * (k*(C-k)))
    return exp(res)
end

function murmur_totality(G, N::Int64=524288, f::Float64=0.2)::Float128
    if f < 0 || f > 1
        throw(ArgumentError("Fraction of Byzantine processes f must be between 0 and 1"))
    end
    C = floor(Int, (1-f)*N)
    fill_log_sum_table(C)
    P::Float128 = compute_p(G, N)
    vals::Array{Float128, 1} = zeros(ceil(Int, C/2))
    for k in 1:ceil(Int, C/2)
        vals[k] = log_sum_body(k, P, C)
    end
    ϵ::Float128 = kahan_summation(vals)
    return check_rounding(ϵ)
end

function run_murmur(G, N::Int64=524290, f::Float64=0.2)
    for n in 1:G
        ϵ = murmur_totality(n)
        println("(", n, ", ", ϵ, ")")
    end
end
