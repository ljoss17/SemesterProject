fp = pwd()*"/Murmur/"

function plot_murmur_G(;N::Int64=102400, f::Float64=0.1, minG::Int64=300, maxG::Int64=1000, step::Int64=1)
    it::Int64 = ((maxG-minG)/step)+1
    x::Array{Int64, 1} = zeros(it)
    y::Array{Float128, 1} = zeros(it)
    f2::Float128 = Float128(f)
    ind = 1
    for g in minG:step:maxG
        x[ind] = g
        y[ind] = murmur_totality(g, N, f2)
        ind = ind+1
    end
    println(x)
    println(y)
    p = plot(
        x,
        y,
        yscale=:log10,
        xlabel="Size of G",
        ylabel="ϵ-totality for Murmur",
        legend=false,
    )
    savefig(p, fp*"murmur_N($N)_f($f).png")
end

function plot_murmur()
    N = 1024
    for n in 10:17
        # Plot for N=1024,2048,...,131'072
        N = 2^n
        plot_murmur_G(N)
    end
end
