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
    d_thr::Int64 = 0
    while left < right
        println("left : $left, right : $right")
        d_thr = floor(Int, (left+right)/2)
        t = contagion_totality(E, E_thr, D, d_thr, R, R_thr, N, f)
        v = contagion_validity(G, E, E_thr, D, d_thr, N, f)
        c = contagion_consistency(E, E_thr, D, d_thr, R, R_thr, N, f)
        ϵ = max(t, v, c)
        println("d_thr : $d_thr, ϵ : $ϵ. t : $t, v : $t, c : $c")
        if ((t+v+c)/3 < ϵ)
            break
        end
        if ϵ == v || ϵ == t
            right = d_thr-1
        else
            left = d_thr+1
        end
    end
    println("D_thr : $d_thr")
    left = 1
    right = d_thr
    while left < right
        println("left : $left, right : $right")
        r_thr = floor(Int, (left+right)/2)
        t = contagion_totality(E, E_thr, D, d_thr, R, r_thr, N, f)
        v = contagion_validity(G, E, E_thr, D, d_thr, N, f)
        c = contagion_consistency(E, E_thr, D, d_thr, R, r_thr, N, f)
        ϵ = max(t, v, c)
        println("r_thr : $r_thr, ϵ : $ϵ. t : $t, v : $t, c : $c")
        if ((t+v+c)/3 < ϵ)
            break
        end
        if ϵ == v || ϵ == t
            right = r_thr-1
        else
            left = r_thr+1
        end
    end
end

function optimize_contagion_R(N::Int64=1024)
    f = 0.1
    G = 18
    E = 150
    E_thr = 103
    D = 150
    d_thr = 95
    R = 150
    left = 1
    right = d_thr
    while left < right
        println("left : $left, right : $right")
        r_thr = floor(Int, (left+right)/2)
        t = contagion_totality(E, E_thr, D, d_thr, R, r_thr, N, f)
        v = contagion_validity(G, E, E_thr, D, d_thr, N, f)
        c = contagion_consistency(E, E_thr, D, d_thr, R, r_thr, N, f)
        ϵ = max(t, v, c)
        println("r_thr : $r_thr, ϵ : $ϵ. t : $t, v : $t, c : $c")
        if ((t+v+c)/3 < ϵ)
            break
        end
        if ϵ == v || ϵ == t
            right = r_thr-1
        else
            left = r_thr+1
        end
    end
end
