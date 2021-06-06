fp = pwd()*"/Sieve/"

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
        title=string(" E_T : ", E_thr),
        xlabel="E",
        ylabel="ϵ-consistency for Sieve",
        legend=false,
    )
    savefig(p, fp*"sieve_consistency_changeE_N($N)_f($f)_E_thr($E_thr)_E_from($minE)_to($maxE)_step($step).png")
end

function plot_sieve_sieve_consistency_E_thr(;N::Int64=1024, f::Float64=0.1, E::Int64=400, minE_thr::Int64=200, maxE_thr::Int64=400, step::Int64=10)
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
        title=string(" E : ", E),
        xlabel="E_thr",
        ylabel="ϵ-consistency for Sieve",
        legend=false,
    )
    savefig(p, fp*"sieve_consistency_changeE-thr_N($N)_f($f)_E($E)_Ethr_from($minE_thr)_to($maxE_thr)_step($step).png")
end

function plot_sieve()
    N = 1024
    E = 300
    E_thr = 100
    G = 100
    step = 1
    for p in 80:10:130
        println("p : $p")
        @time plot_sieve_totalvalidity_E(N=N, G=G, E_thr=p, minE=105, maxE=150, step=1)
        #@time plot_sieve_consistency_E(N=N, E_thr=p, minE=100, maxE=400, step=10)
    end
    for p in 100:10:130
        println("p : $p")
        @time plot_sieve_totalvalidity_E_thr(N=N, G=G, E=p, minE_thr=50, maxE_thr=75, step=1)
        #@time plot_sieve_sieve_consistency_E_thr(N=N, E=p, minE_thr=50, maxE_thr=min(100, p), step=10)
    end
end
