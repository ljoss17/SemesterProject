using Plots

fp = "~/Documents/Julia/SemesterProject/Step1/Sieve/"

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
    x = []
    y = []
    for e in minE:step:maxE
        push!(x, e)
        push!(y, sieve_total_validity(G, e, E_thr))
    end
    p = plot(
        x,
        y,
        title=string("G : ", G, " E_T : ", E_thr),
        xlabel="E",
        ylabel="ϵ-total validity for Sieve",
        legend=false,
    )
    savefig(p, fp*"sieve_totalvalidity_changeE_N($N)_f($f)_G($G)_E_thr($E_thr).png")
end

function plot_sieve_totalvalidity_E_thr(;N::Int64=1024, f::Float64=0.1, G::Int64=300, E::Int64=200, minE_thr::Int64=200, maxE_thr::Int64=400, step::Int64=10)
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
        push!(y, sieve_total_validity(G, E, e_thr))
    end
    p = plot(
        x,
        y,
        title=string("G : ", G, " E : ", E),
        xlabel="E_thr",
        ylabel="ϵ-total validity for Sieve",
        legend=false,
    )
    savefig(p, fp*"sieve_totalvalidity_changeE-thr_N($N)_f($f)_G($G)_E($E).png")
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
    x = []
    y = []
    for e in minE:step:maxE
        push!(x, e)
        push!(y, sieve_consistency(e, E_thr))
    end
    p = plot(
        x,
        y,
        title=string(" E_T : ", E_thr),
        xlabel="E",
        ylabel="ϵ-consistency for Sieve",
        legend=false,
    )
    savefig(p, fp*"sieve_consistency_changeE_N($N)_f($f)_E_thr($E_thr).png")
end

function plot_sieve_sieve_consistency_E_thr(;N::Int64=1024, f::Float64=0.1, E::Int64=200, minE_thr::Int64=200, maxE_thr::Int64=400, step::Int64=10)
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
        push!(y, sieve_consistency(E, e_thr))
    end
    p = plot(
        x,
        y,
        title=string(" E : ", E),
        xlabel="E_thr",
        ylabel="ϵ-consistency for Sieve",
        legend=false,
    )
    savefig(p, fp*"sieve_consistency_changeE-thr_N($N)_f($f)_E($E).png")
end

function plot_sieve()
    N = 1024
    E = 300
    E_thr = 100
    G = 250
    step = 1
    for n in 9:12
        # Plot for N=512,1024,...,4096
        plot_sieve_totalvalidity_E(N=2^n, G=G, E_thr=E_thr, minE=E_thr, maxE=3*E_thr, step=1)
        plot_sieve_consistency_E(N=2^n, E_thr=E_thr, minE=E_thr, maxE=3*E_thr, step=1)
        plot_sieve_totalvalidity_E_thr(N=2^n, G=G, E=E, minE_thr=floor(Int, E/3), maxE_thr=E, step=1)
        plot_sieve_sieve_consistency_E_thr(N=2^n, E=E, minE_thr=floor(Int, E/3), maxE_thr=E, step=1)
    end
end
