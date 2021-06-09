include("contagion.jl")

function compute_epsilon(;N::Int64=1024, f::Float64=0.1, G::Int64=300, E::Int64=300, E_threshold::Int64=100, R::Int64=300, R_threshold::Int64=100, D::Int64=300, D_threshold::Int64=100)
    @time ϵ_t = contagion_totality(E, E_threshold, D, D_threshold, R, R_threshold, N, f)
    @time ϵ_v = contagion_validity(G, E, E_threshold, D, D_threshold, N, f)
    @time ϵ_c = contagion_consistency(E, E_threshold, D, D_threshold, R, R_threshold, N, f)
    return ϵ_t, ϵ_v, ϵ_c
end

function generate_results(n_min::Int64=70, n_max::Int64=130, step::Int64=10)
    println("Start :")
    open("results.txt", "w") do io
        write(io, "For N $n_min to $n_max:\n\n")
    end
    for n in n_min:step:n_max
        println("n : $n")
        ϵ_t, ϵ_v, ϵ_c = compute_epsilon(N=n, G=150, E=150, E_threshold=40, R=150, R_threshold=40, D=150, D_threshold=40)
        ϵ = max(ϵ_t, ϵ_v, ϵ_c)
        open("results.txt", "a") do io
            write(io, "Values for N : $n || ϵ : $ϵ. With ϵ_t : $ϵ_t, ϵ_v : $ϵ_v, ϵ_c : $ϵ_c\n")
        end
    end
end

function generate_results_params(p_min::Int64=150, p_max::Int64=300)
    println("Start :")
    N = 1024
    G = 100
    E = 150
    E_thr = 105
    open("results_params.txt", "w") do io
        write(io, "For N $N:\n\n")
    end
    for p in p_min:10:p_max
        println("p : $p")
        ϵ_t, ϵ_v, ϵ_c = compute_epsilon(N=N, G=G, E=E, E_threshold=E_thr, R=150, R_threshold=p, D=150, D_threshold=70)
        ϵ = max(ϵ_t, ϵ_v, ϵ_c)
        open("results_params.txt", "a") do io
            write(io, "Values for R_thr : $p || ϵ : $ϵ. With ϵ_t : $ϵ_t, ϵ_v : $ϵ_v, ϵ_c : $ϵ_c\n")
        end
        println("R done")
        ϵ_t, ϵ_v, ϵ_c = compute_epsilon(N=N, G=G, E=E, E_threshold=E_thr, R=150, R_threshold=70, D=p, D_threshold=p)
        ϵ = max(ϵ_t, ϵ_v, ϵ_c)
        open("results_params.txt", "a") do io
            write(io, "Values for D_thr : $p || ϵ : $ϵ. With ϵ_t : $ϵ_t, ϵ_v : $ϵ_v, ϵ_c : $ϵ_c\n")
        end
        println("D done")
    end
end
