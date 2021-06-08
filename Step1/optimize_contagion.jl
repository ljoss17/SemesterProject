include("contagion.jl")

function optimize_contagion(N::Int64=1024)
    f = 0.1
    G = 18
    E = 150
    E_thr = 103
    D = 150
    R = 150
    R_thr = 40
    left = R_thr
    right = D
    while left < right
        println("left : $left, right : $right")
        m = floor(Int, (left+right)/2)
        t = contagion_totality(E, E_thr, D, m, R, R_thr, N, f)
        v = contagion_validity(G, E, E_thr, D, m, N, f)
        c = contagion_consistency(E, E_thr, D, m, R, R_thr, N, f)
        ϵ = max(t, v, c)
        if ϵ == v || ϵ == c
            right = m-1
        end
    end
end
