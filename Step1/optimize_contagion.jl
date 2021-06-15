include("contagion.jl")

if Sys.iswindows()
    fp_sieve = pwd()*"\\Parameters\\Sieve\\"
else
    fp_sieve = pwd()*"/Parameters/Sieve/"
end


if Sys.iswindows()
    fp = pwd()*"\\Parameters\\Contagion\\"
else
    fp = pwd()*"/Parameters/Contagion/"
end

function optimize_contagion(N::Int64=1024, E::Int64=150, D::Int64=150, R::Int64=150)
    f = 0.1
    # If Echo is not optimized for the given N and E, compute it.
    if isfile(fp_sieve*"sieve_params_N($N)_E($E).csv")
        file = CSV.File(fp_sieve*"sieve_params_N($N)_E($E).csv")
    else
        sieve_find_threshold(E, E, N)
        file = CSV.File(fp_sieve*"sieve_params_N($N)_E($E).csv")
    end
    G = file.G[1]
    E = file.E[1]
    E_thr = file.E_thr[1]
    leftD = 1
    #rightD = D
    rightD = 99
    d_thr::Int64 = 0
    r_thr::Int64 = 0
    best_ϵ::Float128 = 1.0
    best_v::Float128 = 1.0
    best_c::Float128 = 1.0
    best_t::Float128 = 1.0
    best_e_thr::Int64 = E_thr
    best_d_thr::Int64 = 0
    best_r_thr::Int64 = 0
    ϵ::Float128 = 0.0
    c::Float128 = 0.0
    v::Float128 = 0.0
    t::Float128 = 0.0
    while leftD < rightD
        d_thr = floor(Int, (leftD+rightD)/2)
        println("leftD : $leftD, rightD : $rightD, d_thr : $d_thr")
        leftR = 1
        rightR = d_thr
        v = contagion_validity(G, E, E_thr, D, d_thr, N, f)
        while leftR < rightR
            r_thr = floor(Int, (leftR+rightR)/2)
            println("leftR : $leftR, rightR : $rightR, r_thr : $r_thr")
            t = contagion_totality(E, E_thr, D, d_thr, R, r_thr, N, f)
            c = contagion_consistency(E, E_thr, D, d_thr, R, r_thr, N, f)
            ϵ = max(t, v, c)
            if ϵ < best_ϵ
                best_ϵ = ϵ
                best_v = v
                best_c = c
                best_t = t
                best_d_thr = d_thr
                best_r_thr = r_thr
            end
            # If validity is the limiting factor, changing R_thr won't change ϵ
            # An optimal R_thr for the given D_thr is found.
            if v > t && v > c
                leftR = rightR
                break
            end
            if t > c
                rightR = r_thr-1
            elseif c > t
                leftR = r_thr+1
            else
                break
            end
        end
        open("tmp_res_contagion.txt", "a") do io
            write(io, "d_thr : $d_thr, r_thr : $r_thr, ϵ : $ϵ. t : $t, v : $v, c : $c\n\n")
        end
        println("d_thr : $d_thr, r_thr : $r_thr, ϵ : $ϵ. t : $t, v : $v, c : $c")
        if ϵ == v || ϵ == t
            rightD = d_thr-1
        else
            leftD = d_thr+1
        end
    end
    df = DataFrame(
        ϵ = [best_ϵ],
        ϵ_v = [best_v],
        ϵ_c = [best_c],
        ϵ_t = [best_t],
        N = [N],
        f= [f],
        G = [G],
        E = [E],
        E_thr = [best_e_thr],
        R = [R],
        R_thr = [best_r_thr],
        D = [D],
        D_thr = [best_d_thr]
    )
    CSV.write(fp*"contagion_params_N($N)_E($E)_R($R)_D($D).csv", df)
end
