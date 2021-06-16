include("sieve.jl")

if Sys.iswindows()
    fp = pwd()*"\\Sieve\\"
else
    fp = pwd()*"/Sieve/"
end

function plot_sieve_totalvalidity_E(;N::Int64=1024, f::Float64=0.1, G::Int64=300, E_thr::Int64=200, minE::Int64=200, maxE::Int64=400, step::Int64=10)
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
    it::Int64 = ((maxE-minE)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for e in minE:step:maxE
        x[ind] = e
        y[ind] = sieve_total_validity(G, e, E_thr)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        yscale=:log10,
        title=string("G : ", G, " E_T : ", E_thr),
        xlabel="E",
        ylabel="ϵ-total validity for Sieve",
        legend=false,
    )
    savefig(p, fp*"sieve_totalvalidity_changeE_N($N)_f($f)_G($G)_E_thr($E_thr)_E_from($minE)_to($maxE)_step($step).png")
    return x, y
end

function plot_sieve_totalvalidity_E_thr(;N::Int64=1024, f::Float64=0.1, G::Int64=300, E::Int64=400, minE_thr::Int64=200, maxE_thr::Int64=400, step::Int64=10)
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
    it::Int64 = ((maxE_thr-minE_thr)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for e_thr in minE_thr:step:maxE_thr
        x[ind] = e_thr
        y[ind] = sieve_total_validity(G, E, e_thr)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        yscale=:log10,
        title=string("G : ", G, " E : ", E),
        xlabel="E_thr",
        ylabel="ϵ-total validity for Sieve",
        legend=false,
    )
    savefig(p, fp*"sieve_totalvalidity_changeE-thr_N($N)_f($f)_G($G)_E($E)_Ethr_from($minE_thr)_to($maxE_thr)_step($step).png")
    return x, y
end

function plot_sieve_consistency_E(;N::Int64=1024, f::Float64=0.1, E_thr::Int64=200, minE::Int64=200, maxE::Int64=400, step::Int64=10)
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
    it::Int64 = ((maxE-minE)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for e in minE:step:maxE
        x[ind] = e
        y[ind] = sieve_consistency(e, E_thr)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        yscale=:log10,
        title=string(" E_T : ", E_thr),
        xlabel="E",
        ylabel="ϵ-consistency for Sieve",
        legend=false,
    )
    savefig(p, fp*"sieve_consistency_changeE_N($N)_f($f)_E_thr($E_thr)_E_from($minE)_to($maxE)_step($step).png")
    return x, y
end

function plot_sieve_consistency_E_thr(;N::Int64=1024, f::Float64=0.1, E::Int64=400, minE_thr::Int64=200, maxE_thr::Int64=400, step::Int64=10)
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
    it::Int64 = ((maxE_thr-minE_thr)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    ind = 1
    for e_thr in minE_thr:step:maxE_thr
        x[ind] = e_thr
        y[ind] = sieve_consistency(E, e_thr)
        ind = ind+1
    end
    p = plot(
        x,
        y,
        yscale=:log10,
        title=string(" E : ", E),
        xlabel="E_thr",
        ylabel="ϵ-consistency for Sieve",
        legend=false,
    )
    savefig(p, fp*"sieve_consistency_changeE-thr_N($N)_f($f)_E($E)_Ethr_from($minE_thr)_to($maxE_thr)_step($step).png")
    return x, y
end

function plot_sieve()
    N = 1024
    G = 100
    for p in 300:20:400
        println("p : $p")
        x1, y1 = plot_sieve_totalvalidity_E_thr(N=N, G=G, E=p, minE_thr=p-100, maxE_thr=p-50, step=1)
        x2, y2 = plot_sieve_consistency_E_thr(N=N, E=p, minE_thr=p-100, maxE_thr=p-50, step=1)
        pl = plot(
            x1,
            y1,
            yscale=:log10,
            title=string(" E : ", p),
            xlabel="E_thr",
            ylabel="ϵ-consistency for Sieve",
            labels="Total-validity"
        )
        pl = plot!(
            x2,
            y2,
            yscale=:log10,
            title=string(" E : ", p),
            xlabel="E_thr",
            ylabel="ϵ-consistency for Sieve",
            labels="Consistency"
        )
        display(pl)
        savefig(pl, fp*"sieve_all_p($p)_t.png")
    end
end
