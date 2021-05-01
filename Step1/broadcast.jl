include("contagion.jl")

function compute_epsilon(;N::Int64=1024, f::Float64=0.1, G::Int64=300, E::Int64=300, E_threshold::Int64=100, R::Int64=300, R_threshold::Int64=100, D::Int64=300, D_threshold::Int64=100)
    ϵ_t = contagion_totality(E, E_threshold, D, D_threshold, R, R_threshold, N, f)
    ϵ_v = contagion_validity(G, E, E_threshold, D, D_threshold, N, f)
    ϵ_c = contagion_consistency(E, E_threshold, D, D_threshold, R, R_threshold, N, f)
    return ϵ_t, ϵ_v, ϵ_c
end

function generate_results()
    println("Start :")
    for n in 70:130
        println("n : $n")
        ϵ_t, ϵ_v, ϵ_c = compute_epsilon(N=n, G=15, E=15, E_threshold=5, R=15, R_threshold=5, D=15, D_threshold=5)
        ϵ = max(ϵ_t, ϵ_v, ϵ_c)
        open("results.txt", "a") do io
            write(io, "Value ϵ : $ϵ. With ϵ_t : $ϵ_t, ϵ_v : $ϵ_v, ϵ_c : $ϵ_c\n")
        end
    end
end
