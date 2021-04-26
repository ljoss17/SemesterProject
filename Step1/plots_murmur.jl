using Plots

fp = "~/Documents/Julia/SemesterProject/Step1/Murmur/"

function plot_murmur_G(N::Int64=102400, f::Float64=0.1, minG::Int64=300, maxG::Int64=1000, step::Int64=1)
    x = []
    y = []
    C = floor(Int, (1-f)*N)
    for g in minG:step:maxG
        push!(x, g)
        push!(y, murmur_totality(g, C, N))
    end

    p = plot(
        x,
        y,
        xlabel="Size of G",
        ylabel="ϵ-totality for Murmur",
        legend=false,
    )
    savefig(p, fp*"murmur_N($N)_f($f)")
end

function plot_murmur()
    N = 1024
    for n in 10:17
        # Plot for N=1024,2048,...,131'072
        N = 2^n
        plot_murmur_G(N)
    end
end
