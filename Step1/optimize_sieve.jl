include("sieve.jl")

if Sys.iswindows()
    fp = pwd()*"\\Parameters\\Sieve\\"
else
    fp = pwd()*"/Parameters/Sieve/"
end

function sieve_find_threshold(E::Int64=150, G::Int64=18, N::Int64=1024, f::Float64=0.1)
    # Search for the Echo threshold given a system size, Echo set and Gossip set size.
    left = 1
    right = E
    tmp_ethr::Int64 = 0
    ϵ::Float128 = 1.0
    best_ϵ::Float128 = 1.0
    best_v::Float128 = 1.0
    best_c::Float128 = 1.0
    best_ethr::Int64 = 0
    while left < right
        tmp_ethr = floor(Int, (left+right)/2)
        println("left : $left, right : $right, E_thr : $tmp_ethr")
        v = sieve_total_validity(G, E, tmp_ethr, N, f)
        c = sieve_consistency(E, tmp_ethr, N, f)
        ϵ = max(v, c)
        if ϵ < best_ϵ
            best_ϵ = ϵ
            best_c = c
            best_v = v
            best_ethr = tmp_ethr
        end
        if v < c
            left = tmp_ethr+1
        elseif v > c
            right = tmp_ethr-1
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
        E_thr = [best_ethr]
    )
    CSV.write(fp*"sieve_params_N($N)_E($E).csv", df)
    println("best_ϵ : $best_ϵ, best_g : $best_g, best_ethr : $best_ethr")
    return best_ϵ, best_g, best_ethr
end

function sieve_get_params(bound::Float64=1e-10, N::Int64=1024, f::Float64=0.1)
    # Search for the optimal size of Echo set, Echo threshold and Gossip set given
    # a security bound and system size.
    left = 1
    right = N
    # Murmur totality is computed very quickly so we set G to be very big, since it will be optimized
    # at the same time as the threshold.
    g = N
    ϵ::Float128 = 1.0
    best_ϵ::Float128 = 1.0
    best_e::Int64 = N
    best_g::Int64 = 0
    best_ethr::Int64 = 0
    while left < right
        e = floor(Int, (left+right)/2)
        println("left : $left, right : $right, e : $e")
        tmp_ϵ, tmp_g, tmp_ethr = sieve_find_threshold(e, g, N, f)
        # Keep track of the best ϵ computed, even if desired bound is not reached.
        # If the bound of Sieve is smaller than the desired bound, we can reduce the set E.
        if tmp_ϵ < bound
            # Keep the setup with the smallest E which is under the given bounds.
            if e < best_e
                ϵ = tmp_ϵ
                best_e = e
                best_g = tmp_g
                best_ethr = tmp_ethr
            end
            right = e-1
        else
            left = e+1
        end
    end
    df = DataFrame(
        ϵ = [ϵ],
        N = [N],
        f= [f],
        G = [best_g],
        E = [best_e],
        E_thr = [best_ethr]
    )
    CSV.write(fp*"sieve_params_N($N)_bound($bound).csv", df)
    return ϵ, best_e, best_g, best_ethr
end
