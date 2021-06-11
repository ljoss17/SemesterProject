include("sieve.jl")

if Sys.iswindows()
    fp = pwd()*"\\Parameters\\Sieve\\"
else
    fp = pwd()*"/Parameters/Sieve/"
end

function sieve_find_threshold(E::Int64=150, G::Int64=18, N::Int64=1024)
    f = 0.1
    left = 1
    right = E
    m::Int64 = 0
    ϵ::Float128 = 1.0
    best_ϵ::Float128 = 1.0
    best_v::Float128 = 1.0
    best_c::Float128 = 1.0
    best_e_thr::Int64 = 0
    while left < right
        m = floor(Int, (left+right)/2)
        v = sieve_total_validity(G, E, m, N, f)
        c = sieve_consistency(E, m, N, f)
        ϵ = max(v, c)
        if ϵ < best_ϵ
            best_ϵ = ϵ
            best_c = c
            best_v = v
            best_e_thr = m
        end
        if v < c
            left = m+1
        elseif v > c
            right = m-1
        else
            break
        end
    end
    best_g::Int64 = 0
    for g in 1:E
        t = murmur_totality(g, N, f)
        if t < ϵ*1e-1
            best_g = g
            break
        end
    end
    df = DataFrame(
        ϵ = [best_ϵ],
        ϵ_v = [best_v],
        ϵ_c = [best_c],
        N = [N],
        f= [f],
        G = [best_g],
        E = [E],
        E_thr = [best_e_thr]
    )
    CSV.write(fp*"sieve_params_N($N)_E($E).csv", df)
    return m
end
