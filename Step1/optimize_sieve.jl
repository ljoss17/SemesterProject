include("sieve.jl")

function sieve_find_threshold(E::Int64=150, G::Int64=18, N::Int64=1024)
    f = 0.1
    left = 1
    right = E
    m::Int64 = 0
    ϵ::Float128 = 1.0
    while left < right
        println("left : $left, right : $right")
        m = floor(Int, (left+right)/2)
        v = sieve_total_validity(G, E, m, N, f)
        c = sieve_consistency(E, m, N, f)
        ϵ = max(v, c)
        println("v : $v, c : $c")
        if 3*v < c
            left = m+1
        elseif v > 3*c
            right = m-1
        else
            break
        end
    end
    res_g::Int64 = 0
    for g in 1:E
        t = murmur_totality(g, N, f)
        if t < ϵ*1e-1
            res_g = g
            break
        end
    end
    println("best G : $res_g")
    return m
end
