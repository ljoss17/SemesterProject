using Plots

n = 10

x = []
y = []

for k in 1:n-1
    r = (917^k*4^(n-k) - 916^k*5(n-k))
    push!(x, k)
    push!(y, r)
end

plot(x,y)
