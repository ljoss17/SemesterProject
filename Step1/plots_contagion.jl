include("contagion.jl")

if Sys.iswindows()
    fp = pwd()*"\\Contagion\\"
else
    fp = pwd()*"/Contagion/"
end

function plot_contagion_totality_R(;N::Int64=1024, f::Float64=0.1, E::Int64=300, E_thr::Int64=200, D::Int64=300, D_thr::Int64=200, R_thr::Int64=200, minR::Int64=200, maxR::Int64=400, step::Int64=10)
    if minR > maxR
        println("Error : minimum value of R can't be bigger than max value of R. min R : $minR, max R : $maxR")
        return
    end
    if minR < R_thr
        println("Error : minimum value of R can't be smaller than R_thr. min R : $minR, E_thr : $R_thr")
        return
    end
    if maxR > N
        println("Error : maximum value of R can't be bigger than N. max R : $maxR, N : $N")
        return
    end
    it::Int64 = ((maxR-minR)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for r in minR:step:maxR
        x[ind] = r
        y[ind] = contagion_totality(E, E_thr, D, D_thr, r, R_thr, N, f)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        title=string("E : $E, E_thr : $E_thr, R_thr : $R_thr, D : $D, D_thr : $D_thr"),
        xlabel="R",
        ylabel="ϵ-totality for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_totality_changeR_N($N)_f($f)_E($E)_E_thr($E_thr)_R_thr($R_thr)_D($D)_D_thr($D_thr).png")
    return x, y
end

function plot_contagion_totality_D(;N::Int64=1024, f::Float64=0.1, E::Int64=300, E_thr::Int64=200, D_thr::Int64=200, R::Int64=300, R_thr::Int64=200, minD::Int64=200, maxD::Int64=400, step::Int64=10)
    if minD > maxD
        println("Error : minimum value of D can't be bigger than max value of D. min D : $minD, max E : $maxD")
        return
    end
    if minD < D_thr
        println("Error : minimum value of D can't be smaller than D_thr. min D : $minD, D_thr : $D_thr")
        return
    end
    if maxD > N
        println("Error : maximum value of D can't be bigger than N. max D : $maxD, N : $N")
        return
    end
    it::Int64 = ((maxD-minD)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for d in minD:step:maxD
        x[ind] = d
        y[ind] = contagion_totality(E, E_thr, d, D_thr, R, R_thr, N, f)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        title=string("E : $E, E_thr : $E_thr, R : $R, R_thr : $R_thr, D_thr : $D_thr"),
        xlabel="D",
        ylabel="ϵ-totality for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_totality_changeD_N($N)_f($f)_E($E)_E_thr($E_thr)_R($R)_R_thr($R_thr)_D_thr($D_thr).png")
    return x, y
end

function plot_contagion_totality_R_thr(;N::Int64=1024, f::Float64=0.1, E::Int64=400, E_thr::Int64=200, D::Int64=400, D_thr::Int64=200, R::Int64=400, minR_thr::Int64=200, maxR_thr::Int64=400, step::Int64=10)
    if minR_thr > maxR_thr
        println("Error : minimum value of R can't be bigger than max value of R. min R threshold : $minR_thr, max R threshold : $maxR_thr")
        return
    end
    if maxR_thr > R
        println("Error : maximum value of R threshold can't be bigger than R. max R threshold : $maxR_thr, R : $R")
        return
    end
    if maxR_thr > N
        println("Error : maximum value of R threshold can't be bigger than N. max R threshold : $maxR_thr, N : $N")
        return
    end
    it::Int64 = ((maxR_thr-minR_thr)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for r_thr in minR_thr:step:maxR_thr
        x[ind] = r_thr
        y[ind] = contagion_totality(E, E_thr, D, D_thr, R, r_thr, N, f)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        title=string("E : $E, E_thr : $E_thr, R : $R, D : $D, D_thr : $D_thr"),
        xlabel="R_thr",
        ylabel="ϵ-totality for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_totality_changeR-thr_N($N)_f($f)_E($E)_E_thr($E_thr)_R($R)_D($D)_D_thr($D_thr).png")
    return x, y
end

function plot_contagion_totality_D_thr(;N::Int64=1024, f::Float64=0.1, E::Int64=400, E_thr::Int64=200, D::Int64=400, R::Int64=400, R_thr::Int64=200, minD_thr::Int64=200, maxD_thr::Int64=400, step::Int64=10)
    if minD_thr > maxD_thr
        println("Error : minimum value of D can't be bigger than max value of D. min D threshold : $minD_thr, max D threshold : $maxD_thr")
        return
    end
    if maxD_thr > D
        println("Error : maximum value of D threshold can't be bigger than D. max D threshold : $maxD_thr, D : $D")
        return
    end
    if maxD_thr > N
        println("Error : maximum value of D threshold can't be bigger than N. max D threshold : $maxD_thr, N : $N")
        return
    end
    it::Int64 = ((maxD_thr-minD_thr)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for d_thr in minD_thr:step:maxD_thr
        x[ind] = d_thr
        y[ind] = contagion_totality(E, E_thr, D, d_thr, R, R_thr, N, f)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        title=string("E : $E, E_thr : $E_thr, R : $R, R_thr : $R_thr, D : $D"),
        xlabel="D_thr",
        ylabel="ϵ-totality for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_totality_changeD-thr_N($N)_f($f)_E($E)_E_thr($E_thr)_R($R)_R_thr($R_thr)_D($D).png")
    return x, y
end

function plot_contagion_validity_D(;N::Int64=1024, f::Float64=0.1, G::Int64=300, E::Int64=300, E_thr::Int64=200, D_thr::Int64=200, minD::Int64=200, maxD::Int64=400, step::Int64=10)
    if minD > maxD
        println("Error : minimum value of D can't be bigger than max value of D. min D : $minD, max E : $maxD")
        return
    end
    if minD < D_thr
        println("Error : minimum value of D can't be smaller than D_thr. min D : $minD, D_thr : $D_thr")
        return
    end
    if maxD > N
        println("Error : maximum value of D can't be bigger than N. max D : $maxD, N : $N")
        return
    end
    it::Int64 = ((maxD-minD)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for d in minD:step:maxD
        x[ind] = d
        y[ind] = contagion_validity(G, E, E_thr, d, D_thr, N, f)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        title=string("G : $G, E : $E, E_thr : $E_thr, D_thr : $D_thr"),
        xlabel="D",
        ylabel="ϵ-validity for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_validity_changeD_N($N)_f($f)_G($G)_E($E)_E_thr($E_thr)_D_thr($D_thr).png")
    return x, y
end

function plot_contagion_validity_D_thr(;N::Int64=1024, f::Float64=0.1, G::Int64=300, E::Int64=400, E_thr::Int64=200, D::Int64=400, minD_thr::Int64=200, maxD_thr::Int64=400, step::Int64=10)
    if minD_thr > maxD_thr
        println("Error : minimum value of D can't be bigger than max value of D. min D threshold : $minD_thr, max D threshold : $maxD_thr")
        return
    end
    if maxD_thr > D
        println("Error : maximum value of D threshold can't be bigger than D. max D threshold : $maxD_thr, D : $D")
        return
    end
    if maxD_thr > N
        println("Error : maximum value of D threshold can't be bigger than N. max D threshold : $maxD_thr, N : $N")
        return
    end
    it::Int64 = ((maxD_thr-minD_thr)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for d_thr in minD_thr:step:maxD_thr
        x[ind] = d_thr
        y[ind] = contagion_validity(G, E, E_thr, D, d_thr, N, f)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        title=string("G : $G, E : $E, E_thr : $E_thr, D : $D"),
        xlabel="D_thr",
        ylabel="ϵ-validity for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_validity_changeD-thr_N($N)_f($f)_G($G)_E($E)_E_thr($E_thr)_D($D).png")
    return x, y
end

function plot_contagion_consistency_R(;N::Int64=1024, f::Float64=0.1, E::Int64=300, E_thr::Int64=200, D::Int64=300, D_thr::Int64=200, R_thr::Int64=200, minR::Int64=200, maxR::Int64=400, step::Int64=10)
    if minR > maxR
        println("Error : minimum value of R can't be bigger than max value of R. min R : $minR, max R : $maxR")
        return
    end
    if minR < R_thr
        println("Error : minimum value of R can't be smaller than R_thr. min R : $minR, E_thr : $R_thr")
        return
    end
    if maxR > N
        println("Error : maximum value of R can't be bigger than N. max R : $maxR, N : $N")
        return
    end
    it::Int64 = ((maxR-minR)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for r in minR:step:maxR
        x[ind] = r
        y[ind] = contagion_consistency(E, E_thr, D, D_thr, r, R_thr, N, f)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        yscale=:log10,
        title=string("E : $E, E_thr : $E_thr, R_thr : $R_thr, D : $D, D_thr : $D_thr"),
        xlabel="R",
        ylabel="ϵ-consistency for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_consistency_changeR_N($N)_f($f)_E($E)_E_thr($E_thr)_R_thr($R_thr)_D($D)_D_thr($D_thr).png")
    return x, y
end

function plot_contagion_consistency_D(;N::Int64=1024, f::Float64=0.1, E::Int64=300, E_thr::Int64=200, D_thr::Int64=200, R::Int64=300, R_thr::Int64=200, minD::Int64=200, maxD::Int64=400, step::Int64=10)
    if minD > maxD
        println("Error : minimum value of D can't be bigger than max value of D. min D : $minD, max E : $maxD")
        return
    end
    if minD < D_thr
        println("Error : minimum value of D can't be smaller than D_thr. min D : $minD, D_thr : $D_thr")
        return
    end
    if maxD > N
        println("Error : maximum value of D can't be bigger than N. max D : $maxD, N : $N")
        return
    end
    it::Int64 = ((maxD-minD)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for d in minD:step:maxD
        x[ind] = d
        y[ind] = contagion_consistency(E, E_thr, d, D_thr, R, R_thr, N, f)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        title=string("E : $E, E_thr : $E_thr, R : $R, R_thr : $R_thr, D_thr : $D_thr"),
        xlabel="D",
        ylabel="ϵ-consistency for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_consistency_changeD_N($N)_f($f)_E($E)_E_thr($E_thr)_R($R)_R_thr($R_thr)_D_thr($D_thr).png")
    return x, y
end

function plot_contagion_consistency_R_thr(;N::Int64=1024, f::Float64=0.1, E::Int64=400, E_thr::Int64=200, D::Int64=400, D_thr::Int64=200, R::Int64=400, minR_thr::Int64=200, maxR_thr::Int64=400, step::Int64=10)
    if minR_thr > maxR_thr
        println("Error : minimum value of R can't be bigger than max value of R. min R threshold : $minR_thr, max R threshold : $maxR_thr")
        return
    end
    if maxR_thr > R
        println("Error : maximum value of R threshold can't be bigger than R. max R threshold : $maxR_thr, R : $R")
        return
    end
    if maxR_thr > N
        println("Error : maximum value of R threshold can't be bigger than N. max R threshold : $maxR_thr, N : $N")
        return
    end
    it::Int64 = ((maxR_thr-minR_thr)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for r_thr in minR_thr:step:maxR_thr
        x[ind] = r_thr
        y[ind] = contagion_consistency(E, E_thr, D, D_thr, R, r_thr, N, f)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        title=string("E : $E, E_thr : $E_thr, R : $R, D : $D, D_thr : $D_thr"),
        xlabel="R_thr",
        ylabel="ϵ-consistency for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_consistency_changeR-thr_N($N)_f($f)_E($E)_E_thr($E_thr)_R($R)_D($D)_D_thr($D_thr).png")
    return x, y
end

function plot_contagion_consistency_D_thr(;N::Int64=1024, f::Float64=0.1, E::Int64=400, E_thr::Int64=200, D::Int64=400, R::Int64=400, R_thr::Int64=200, minD_thr::Int64=200, maxD_thr::Int64=400, step::Int64=10)
    if minD_thr > maxD_thr
        println("Error : minimum value of D can't be bigger than max value of D. min D threshold : $minD_thr, max D threshold : $maxD_thr")
        return
    end
    if maxD_thr > D
        println("Error : maximum value of D threshold can't be bigger than D. max D threshold : $maxD_thr, D : $D")
        return
    end
    if maxD_thr > N
        println("Error : maximum value of D threshold can't be bigger than N. max D threshold : $maxD_thr, N : $N")
        return
    end
    it::Int64 = ((maxD_thr-minD_thr)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for d_thr in minD_thr:step:maxD_thr
        x[ind] = d_thr
        y[ind] = contagion_consistency(E, E_thr, D, d_thr, R, R_thr, N, f)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        title=string("E : $E, E_thr : $E_thr, R : $R, R_thr : $R_thr, D : $D"),
        xlabel="D_thr",
        ylabel="ϵ-consistency for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_consistency_changeD-thr_N($N)_f($f)_E($E)_E_thr($E_thr)_R($R)_R_thr($R_thr)_D($D).png")
    return x, y
end

function plot_gamma_distribution(;N::Int64=60, R::Int64=30, R_t::Int64=12, f::Float64=0.1, minX::Int64=115, maxX::Int64=128)
    C = floor(Int, (1-f)*N)
    S = N-C
    γ_end = N
    γs = threshold_contagion(N, R, 1, 1, S, R_t, γ_end)
    if length(γs) < maxX+1
    end
    println(γs)
    p = plot(
        [minX:maxX],
        γs[minX+1:maxX+1]
    )
    display(p)
    savefig(p, fp*"/gammas/N_($N)")
end

function plot_consistency_gammas(;r::Int64=35, r_t::Int64=12, minN::Int64=100, maxN::Int64=300)
    minC = floor(Int, (1-0.1)*minN)
    maxN = 300
    for n in minN:10:maxN
        plot_gamma_distribution(N=n, R=r, R_t=r_t, minX=minC, maxX=maxN)
    end
end

function plot_totality_gammas(;N::Int64=1024, r::Int64=200, r_t::Int64=60)
    C = floor(Int, 0.9*N)
    dists = contagion_threshold(C, r, 0.9, C, 1, r_t, C)
    k = 0
    for y in dists
        k = k+1
        p = plot(
            [1:C+1],
            y,
            ylims=[0,1],
            title=string("N : $C, S : 1, R : $r, K : $k, l : 0.9, R_thr : $r_t"),
            xlabel="γ+",
            ylabel="p[γ(N, S, R, K, l, R_thr)=γ+]",
            legend=false
        )
        savefig(p, fp*"gammas_totality_K_($k)")
    end
end

function plot_contagion_allD()
    for p in 40:10:80
        x1, y1 = plot_contagion_consistency_D_thr(E=250, E_thr=197, D=200, R=200, R_thr=p, minD_thr=p, maxD_thr=p+120, step=1)
        x2, y2 = plot_contagion_totality_D_thr(E=250, E_thr=197, D=200, R=200, R_thr=p, minD_thr=p, maxD_thr=p+120, step=1)
        x3, y3 = plot_contagion_validity_D_thr(G=11, E=250, E_thr=197, D=200, minD_thr=p, maxD_thr=p+120, step=1)
        plot(
            x1,
            y1,
            yscale=:log10,
            title=string(" D : $D, R : $R"),
            xlabel="E_thr",
            ylabel="ϵ-consistency for Sieve",
            labels="Consistency"
        )
        plot!(
            x2,
            y2,
            yscale=:log10,
            title=string(" D : $D, R : $R"),
            xlabel="E_thr",
            ylabel="ϵ-consistency for Sieve",
            labels="Totality"
        )

        plot!(
            x3,
            y3,
            yscale=:log10,
            title=string(" D : $D, R : $R"),
            xlabel="E_thr",
            ylabel="ϵ-consistency for Sieve",
            labels="Validity"
        )
        savefig(fp*"sieve_all_p($p).png")
    end
end
