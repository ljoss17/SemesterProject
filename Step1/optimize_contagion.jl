include("contagion.jl")
include("optimize_sieve.jl")

if Sys.iswindows()
    fp_c = pwd()*"\\Parameters\\Contagion\\"
else
    fp_c = pwd()*"/Parameters/Contagion/"
end

function optimize_contagion_thresholds(N::Int64=1024, f::Float64=0.1, G::Int64=100, E::Int64=150, E_thr::Int64=100, D::Int64=150, R::Int64=150)
    # Search for the optimal Ready and Deliver threshold, given all the other parameters.
    leftD::Int64 = 1
    rightD::Int64 = D
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
    next_t_r::Float128 = 1.0
    next_best_t_d::Float128 = 1.0
    while leftD < rightD
        d_thr = floor(Int, (leftD+rightD)/2)
        println("leftD : $leftD, rightD : $rightD, d_thr : $d_thr")
        leftR = 1
        # The requirements are R_thr/R < D_thr/D and  R_thr < R.
        rightR = min(floor(Int, (d_thr*R)/D), R)
        v = contagion_validity(G, E, E_thr, D, d_thr, N, f)
        while leftR < rightR
            r_thr = floor(Int, (leftR+rightR)/2)
            println("leftR : $leftR, rightR : $rightR, r_thr : $r_thr")
            t = contagion_totality(E, E_thr, D, d_thr, R, r_thr, N, f)
            # If R_thr can be incremented by 1, compute to observe if t is increasing or decreasing
            # Else force the decreasing of R_thr.
            if r_thr < R && r_thr < floor(Int, (d_thr*R)/D)
                next_t_r = contagion_totality(E, E_thr, D, d_thr, R, r_thr+1, N, f)
            else
                next_t_r = 1.1
            end
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
            # An optimal R_thr for the given D_thr has been found.
            if v > t && v > c
                leftR = rightR
                break
            end
            # If totality is higher, observe if totality is increasing or decreasing.
            if t > c
                # If totality is increasing, decrease threshold
                # If totality is decreasing, increase threshold
                if t > next_t_r
                    leftR = r_thr+1
                else
                    rightR = r_thr-1
                end
            # If consistency is higher, increase the threshold
            else
                leftR = r_thr+1
            end
        end
        open("tmp_res_contagion.txt", "a") do io
            write(io, "d_thr : $d_thr, r_thr : $r_thr, ϵ : $ϵ. t : $t, v : $v, c : $c\n\n")
        end
        # If totality is the highest bound, check if it's increaseing or decreasing
        if best_ϵ == best_t
            # If D_thr can be incremented by 1, check if totality is increasing or decreasing.
            # Else decrease set the next_best_t_d to decrease D_thr.
            if d_thr < D
                next_best_t_d = contagion_totality(E, E_thr, D, d_thr+1, R, r_thr, N, f)
            else
                next_best_t_d = 1.1
            end
            if best_t < next_best_t_d
                rightD = d_thr-1
            else
                leftD = d_thr+1
            end
        # Else check if consistency is higher or lower than validity
        else
            if best_c < best_v
                rightD = d_thr-1
            elseif best_c > best_v
                leftD = d_thr+1
            else
                leftD = rightD
                break
            end
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
    CSV.write(fp_c*"contagion_params_N($N)_G($G)_Ethr($E_thr)_E($E)_R($R)_D($D).csv", df)
    println("Best for thresholds : $best_ϵ")
    return best_ϵ, best_r_thr, best_d_thr
end

function optimize_contagion(N::Int64=1024, f::Float64=0.1, bound::Float64=1e-10)
    # Search for all the optimal parameters for a given system size and security bound.
    leftD::Int64 = 1
    rightD::Int64 = N
    ϵ, e, g, ethr = sieve_get_params(bound*1e-1, N, f)
    tmp_d::Int64 = 0
    tmp_r::Int64 = 0
    leftR::Int64 = 0
    rightR::Int64 = 0
    ϵ_d::Float128 = 1.0
    ϵ_r::Float128 = 1.0
    best_ϵ::Float128 = 1.0
    best_d::Int64 = 0
    intermerdiary_r::Int64 = 0
    best_r::Int64 = 0
    best_dthr::Int64 = 0
    intermerdiary_rthr::Int64 = 0
    best_rthr::Int64 = 0
    intermerdiary_dthr::Int64 = 0
    tmp_avg::Int64 = N
    while leftD < rightD
        tmp_d = floor(Int, (leftD+rightD)/2)
        println("leftD : $leftD, rightD : $rightD, D : $tmp_d")
        leftR = 1
        rightR = N
        ϵ_d = 1.0
        while leftR < rightR
            tmp_r = floor(Int, (leftR+rightR)/2)
            println("leftR : $leftR, rightR : $rightR, R : $tmp_r")
            ϵ_r, tmp_rthr, tmp_dthr = optimize_contagion_thresholds(N, f, g, e, ethr, tmp_d, tmp_r)
            if ϵ_r < bound
                if ϵ_r < ϵ_d
                    ϵ_d = ϵ_r
                    intermerdiary_r = tmp_r
                    intermerdiary_rthr = tmp_rthr
                    intermerdiary_dthr = tmp_dthr
                end
                rightR = tmp_r-1
            else
                leftR = tmp_r+1
            end
        end
        if ϵ_d < bound
            act_avg = ceil(Int, (tmp_d+intermerdiary_r)/2)
            if act_avg < tmp_avg
                tmp_avg = act_avg
                best_ϵ = ϵ_d
                best_d = tmp_d
                best_dthr = intermerdiary_dthr
                best_r = intermerdiary_r
                best_rthr = intermerdiary_rthr
            end
            rightD = tmp_d-1
        else
            leftD = tmp_d+1
        end
    end
    avg::Int64 = ceil(Int, (g+e+best_r+best_d)/4)
    df = DataFrame(
        ϵ = [best_ϵ],
        N = [N],
        f= [f],
        G = [g],
        E = [e],
        E_thr = [ethr],
        R = [best_r],
        R_thr = [best_rthr],
        D = [best_d],
        D_thr = [best_dthr],
        avg = [avg]
    )
    CSV.write(fp_c*"contagion_params_N($N)_bound($bound).csv", df)
end
