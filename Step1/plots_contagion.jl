include("contagion.jl")

fp = pwd()*"\\Contagion\\"

function plot_contagion_totality_E(;N::Int64=1024, f::Float64=0.1, E_thr::Int64=200, D::Int64=300, D_thr::Int64=200, R::Int64=300, R_thr::Int64=200, minE::Int64=200, maxE::Int64=400, step::Int64=10)
    if minE > maxE
        println("Error : minimum value of E can't be bigger than max value of E. min E : $minE, max E : $maxE")
        return
    end
    if minE < E_thr
        println("Error : minimum value of E can't be smaller than E_thr. min E : $minE, E_thr : $E_thr")
        return
    end
    if maxE > N
        println("Error : maximum value of E can't be bigger than N. max E : $maxE, N : $N")
        return
    end
    x = []
    y = []
    for e in minE:step:maxE
        push!(x, e)
        push!(y, contagion_totality(e, E_thr, D, D_thr, R, R_thr, N, f))
    end
    p = plot(
        x,
        y,
        title=string("E_thr : $E_thr, R : $R, R_thr : $R_thr, D : $D, D_thr : $D_thr"),
        xlabel="E",
        ylabel="ϵ-totality for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_totality_changeE_N($N)_f($f)_E_thr($E_thr)_R($R)_R_thr($R_thr)_D($D)_D_thr($D_thr).png")
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
    x = []
    y = []
    for r in minR:step:maxR
        push!(x, r)
        push!(y, contagion_totality(E, E_thr, D, D_thr, r, R_thr, N, f))
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
    x = []
    y = []
    for d in minD:step:maxD
        push!(x, d)
        push!(y, contagion_totality(E, E_thr, d, D_thr, R, R_thr, N, f))
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
end

function plot_contagion_totality_E_thr(;N::Int64=1024, f::Float64=0.1, E::Int64=400, D::Int64=400, D_thr::Int64=200, R::Int64=400, R_thr::Int64=200, minE_thr::Int64=200, maxE_thr::Int64=400, step::Int64=10)
    if minE_thr > maxE_thr
        println("Error : minimum value of E can't be bigger than max value of E. min E threshold : $minE_thr, max E threshold : $maxE_thr")
        return
    end
    if maxE_thr > E
        println("Error : maximum value of E threshold can't be bigger than E. max E threshold : $maxE_thr, E : $E")
        return
    end
    if maxE_thr > N
        println("Error : maximum value of E threshold can't be bigger than N. max E threshold : $maxE_thr, N : $N")
        return
    end
    x = []
    y = []
    for e_thr in minE_thr:step:maxE_thr
        push!(x, e_thr)
        push!(y, contagion_totality(E, e_thr, D, D_thr, R, R_thr, N, f))
    end
    p = plot(
        x,
        y,
        title=string("E : $E, R : $R, R_thr : $R_thr, D : $D, D_thr : $D_thr"),
        xlabel="E_thr",
        ylabel="ϵ-totality for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_totality_changeE-thr_N($N)_f($f)_E($E)_R($R)_R_thr($R_thr)_D($D)_D_thr($D_thr).png")
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
    x = []
    y = []
    for r_thr in minR_thr:step:maxR_thr
        push!(x, r_thr)
        push!(y, contagion_totality(E, E_thr, D, D_thr, R, r_thr, N, f))
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
    x = []
    y = []
    for d_thr in minD_thr:step:maxD_thr
        push!(x, d_thr)
        push!(y, contagion_totality(E, E_thr, D, d_thr, R, R_thr, N, f))
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
end

function plot_contagion_validity_G(;N::Int64=1024, f::Float64=0.1, E::Int64=300, E_thr::Int64=200, D::Int64=300, D_thr::Int64=200, minG::Int64=200, maxG::Int64=400, step::Int64=10)
    if minG > maxG
        println("Error : minimum value of G can't be bigger than max value of G. min G : $minG, max G : $maxG")
        return
    end
    if maxG > N
        println("Error : maximum value of G can't be bigger than N. max G : $maxG, N : $N")
        return
    end
    x = []
    y = []
    for g in minG:step:maxG
        push!(x, g)
        push!(y, contagion_validity(g, E, E_thr, D, D_thr, N, f))
    end
    p = plot(
        x,
        y,
        yscale=:log10,
        title=string("E : $E, E_thr : $E_thr, D : $D, D_thr : $D_thr"),
        xlabel="G",
        ylabel="ϵ-validity for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_validity_changeG_N($N)_f($f)_E($E)_E_thr($E_thr)_D($D)_D_thr($D_thr).png")
end

function plot_contagion_validity_E(;N::Int64=1024, f::Float64=0.1, G::Int64=300, E_thr::Int64=200, D::Int64=300, D_thr::Int64=200, minE::Int64=200, maxE::Int64=400, step::Int64=10)
    if minE > maxE
        println("Error : minimum value of E can't be bigger than max value of E. min E : $minE, max E : $maxE")
        return
    end
    if minE < E_thr
        println("Error : minimum value of E can't be smaller than E_thr. min E : $minE, E_thr : $E_thr")
        return
    end
    if maxE > N
        println("Error : maximum value of E can't be bigger than N. max E : $maxE, N : $N")
        return
    end
    x = []
    y = []
    for e in minE:step:maxE
        push!(x, e)
        push!(y, contagion_validity(G, e, E_thr, D, D_thr, N, f))
    end
    p = plot(
        x,
        y,
        yscale=:log10,
        title=string("G : $G, E_thr : $E_thr, D : $D, D_thr : $D_thr"),
        xlabel="E",
        ylabel="ϵ-validity for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_validity_changeE_N($N)_f($f)_G($G)_E_thr($E_thr)_D($D)_D_thr($D_thr).png")
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
    x = []
    y = []
    for d in minD:step:maxD
        push!(x, d)
        push!(y, contagion_validity(G, E, E_thr, d, D_thr, N, f))
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
end

function plot_contagion_validity_E_thr(;N::Int64=1024, f::Float64=0.1, G::Int64=300, E::Int64=400, D::Int64=400, D_thr::Int64=200, minE_thr::Int64=200, maxE_thr::Int64=400, step::Int64=10)
    if minE_thr > maxE_thr
        println("Error : minimum value of E can't be bigger than max value of E. min E threshold : $minE_thr, max E threshold : $maxE_thr")
        return
    end
    if maxE_thr > E
        println("Error : maximum value of E threshold can't be bigger than E. max E threshold : $maxE_thr, E : $E")
        return
    end
    if maxE_thr > N
        println("Error : maximum value of E threshold can't be bigger than N. max E threshold : $maxE_thr, N : $N")
        return
    end
    x = []
    y = []
    for e_thr in minE_thr:step:maxE_thr
        push!(x, e_thr)
        push!(y, contagion_validity(G, E, e_thr, D, D_thr, N, f))
    end
    p = plot(
        x,
        y,
        title=string("G : $G, E : $E, D : $D, D_thr : $D_thr"),
        xlabel="E_thr",
        ylabel="ϵ-validity for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_validity_changeE-thr_N($N)_f($f)_G($G)_E($E)_D($D)_D_thr($D_thr).png")
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
    x = []
    y = []
    for d_thr in minD_thr:step:maxD_thr
        push!(x, d_thr)
        push!(y, contagion_validity(G, E, E_thr, D, d_thr, N, f))
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
end


function plot_contagion_consistency_E(;N::Int64=1024, f::Float64=0.1, E_thr::Int64=200, D::Int64=300, D_thr::Int64=200, R::Int64=300, R_thr::Int64=200, minE::Int64=200, maxE::Int64=400, step::Int64=10)
    if minE > maxE
        println("Error : minimum value of E can't be bigger than max value of E. min E : $minE, max E : $maxE")
        return
    end
    if minE < E_thr
        println("Error : minimum value of E can't be smaller than E_thr. min E : $minE, E_thr : $E_thr")
        return
    end
    if maxE > N
        println("Error : maximum value of E can't be bigger than N. max E : $maxE, N : $N")
        return
    end
    x = []
    y = []
    for e in minE:step:maxE
        push!(x, e)
        push!(y, contagion_consistency(e, E_thr, D, D_thr, R, R_thr, N, f))
    end
    p = plot(
        x,
        y,
        title=string("E_thr : $E_thr, R : $R, R_thr : $R_thr, D : $D, D_thr : $D_thr"),
        xlabel="E",
        ylabel="ϵ-consistency for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_consistency_changeE_N($N)_f($f)_E_thr($E_thr)_R($R)_R_thr($R_thr)_D($D)_D_thr($D_thr).png")
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
    x = []
    y = []
    for r in minR:step:maxR
        push!(x, r)
        push!(y, contagion_consistency(E, E_thr, D, D_thr, r, R_thr, N, f))
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
    x = []
    y = []
    for d in minD:step:maxD
        push!(x, d)
        push!(y, contagion_consistency(E, E_thr, d, D_thr, R, R_thr, N, f))
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
end

function plot_contagion_consistency_E_thr(;N::Int64=1024, f::Float64=0.1, E::Int64=400, D::Int64=400, D_thr::Int64=200, R::Int64=400, R_thr::Int64=200, minE_thr::Int64=200, maxE_thr::Int64=400, step::Int64=10)
    if minE_thr > maxE_thr
        println("Error : minimum value of E can't be bigger than max value of E. min E threshold : $minE_thr, max E threshold : $maxE_thr")
        return
    end
    if maxE_thr > E
        println("Error : maximum value of E threshold can't be bigger than E. max E threshold : $maxE_thr, E : $E")
        return
    end
    if maxE_thr > N
        println("Error : maximum value of E threshold can't be bigger than N. max E threshold : $maxE_thr, N : $N")
        return
    end
    x = []
    y = []
    for e_thr in minE_thr:step:maxE_thr
        push!(x, e_thr)
        push!(y, contagion_consistency(E, e_thr, D, D_thr, R, R_thr, N, f))
    end
    p = plot(
        x,
        y,
        title=string("E : $E, R : $R, R_thr : $R_thr, D : $D, D_thr : $D_thr"),
        xlabel="E_thr",
        ylabel="ϵ-consistency for Contagion",
        legend=false,
    )
    savefig(p, fp*"contagion_consistency_changeE-thr_N($N)_f($f)_E($E)_R($R)_R_thr($R_thr)_D($D)_D_thr($D_thr).png")
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
    x = []
    y = []
    for r_thr in minR_thr:step:maxR_thr
        push!(x, r_thr)
        push!(y, contagion_consistency(E, E_thr, D, D_thr, R, r_thr, N, f))
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
    x = []
    y = []
    for d_thr in minD_thr:step:maxD_thr
        println("d_thr : $d_thr")
        push!(x, d_thr)
        push!(y, contagion_consistency(E, E_thr, D, d_thr, R, R_thr, N, f))
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
end

# Testing plots

function plot_mu()
    x = []
    y = []
    N = 332
    D = 80
    D_t = 55
    C = floor(Int, 0.9*N)
    for u in N-C:N
        v = (1-(1-specific_enough_ready(D, D_t, u, N))^C)
        push!(x, u)
        push!(y, v)
    end
    display(
        plot(
            x,
            y,
            #yscale = :log10,
            title=string("N : ", N, " D : ", D, " D_thr : ", D_t),
            xlabel="gamma",
            ylabel="μ",
            legend=false,
        )
    )
end

function plot_gamma_consistency(;N::Int64=1024, f::Float64=0.1, R::Int64=50, R_t::Int64=25, γ_min::Int64=103, γ_max::Int64=1024, pre::String="")
    x = []
    y = []
    C = floor(Int, (1-f)*N)

    for u in γ_min:γ_max
        v = compute_gamma(N, R, 1, 1, N-C, R_t, u)
        push!(x, u)
        push!(y, v)
    end

    p = plot(
        x,
        y,
        title=string("N : ", N, " R : ", R, " R_thr : ", R_t),
        xlabel="gamma",
        ylabel="p[γ]",
        legend=false,
    )
    savefig(p, fp*pre*"gamma_consistency_N($N)_R($R)_R_t($R_t)")
end

function plot_gamma_totality(;N::Int64=1024, f::Float64=0.1, K::Int64=921, R::Int64=50, R_t::Int64=25, γ_min::Int64=0, γ_max::Int64=921, pre::String="")
    x = []
    y = []
    C = floor(Int, 0.9*N)

    for u in γ_min:γ_max
        v = compute_gamma(C, R, 1-f, K, 1, R_t, u)
        push!(x, u)
        push!(y, v)
    end
    p = plot(
        x,
        y,
        title=string("N : ", N, " R : ", R, " R_thr : ", R_t),
        xlabel="gamma",
        ylabel="p[γ]",
        legend=false,
    )
    savefig(p, fp*pre*"gamma_totality_N($N)_K($K)_R($R)_R_t($R_t)")
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
            title=string("Gamma distribtuion, round : $k"),
            xlabel="gamma",
            ylabel="p[γ]",
            legend=false
        )
        savefig(p, fp*"gammas_totality_K_($k)")
    end

end

function run_gammas()
    R_thr = 100
    for r in 100:50:300
        println("r : $r")
        plot_gamma_consistency(R=r, R_t=R_thr)
        plot_gamma_totality(R=r, R_t=R_thr, K=1)
        plot_gamma_totality(R=r, R_t=R_thr, K=200)
        plot_gamma_totality(R=r, R_t=R_thr, K=500)
        plot_gamma_totality(R=r, R_t=R_thr, K=921)
    end
end

function run_gammas_focused()
    R_thr = 100
    for r in 100:50:300
        println("rf : $r")
        plot_gamma_consistency(R=r, R_t=R_thr, γ_min=0, γ_max=150, pre="focus_")
        plot_gamma_totality(R=r, R_t=R_thr, K=1, γ_min=0, γ_max=30, pre="focus_")
        plot_gamma_totality(R=r, R_t=R_thr, K=200, γ_min=180, γ_max=220, pre="focus_")
        plot_gamma_totality(R=r, R_t=R_thr, K=400, γ_min=380, γ_max=420, pre="focus_")
        plot_gamma_totality(R=r, R_t=R_thr, K=600, γ_min=580, γ_max=620, pre="focus_")
        plot_gamma_totality(R=r, R_t=R_thr, K=921, γ_min=100, γ_max=150, pre="focus_")
    end
end

function run_all()
    run_gammas()
    run_gammas_focused()
end
