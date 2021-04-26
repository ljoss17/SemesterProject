using Plots
include("contagion.jl")

fp = "~/Documents/Julia/SemesterProject/Step1/"

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

function plot_gamma(N::Int64=355, R::Int64=50, R_t::Int64=25, add::String="test")
    x = []
    y = []
    C = floor(Int, 0.9*N)

    for u in N-C:N
        v = compute_gamma(N, R, 1, 1, N-C, R_t, u, 0, 0)
        push!(x, u)
        push!(y, v)
    end
    p = plot(
        x,
        y,
        #yscale = :log10,
        title=string("N : ", N, " D : ", R, " R_thr : ", R_t),
        xlabel="gamma",
        ylabel="p[γ]",
        legend=false,
    )
    savefig(p, fp*add*"gamma_N($N)_R($R)_R_t($R_t)")
end

function run_gammas()
    for v in 200:50:500
        println("v : $v")
        n = 400
        r = 100
        r_t = 66
        plot_gamma(v, r, r_t, "N/")
        plot_gamma(n, v-100, r_t, "R/")
        plot_gamma(n, r, floor(Int,(v/4)-25), "R_thr/")
    end
end
