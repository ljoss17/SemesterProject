using Plots

x = []
y = []
p = []
N = 102400
f = 0.2
C = floor(Int, (1-f)*N)
for g in 100:10:300
    push!(x, g)
    push!(y, murmur_totality(g, C, N))
    push!(p, compute_p(g, N))
end

p1 = plot(
    x,
    y,
    yscale = :log10,
    xlabel="Size of G",
    ylabel="ϵ-totality for Murmur",
    legend=false,
)
p2 = plot(
    x,
    p,
    yscale = :log10,
    xlabel="Size of G",
    ylabel="Erdős–Rényi edge probability p",
    legend=false,
)
plot(p1, p2)
