include("contagion.jl")

function compute_epsilon(N::Int64=1024, f::Float64=0.1, G::Int64=300, E::Int64=300, E_threshold::Int64=100, R::Int64=300, R_threshold::Int64=100, D::Int64=300, D_threshold::Int64=100)
    ϵ_t = contagion_totality(E, E_threshold, D, D_threshold, R, R_threshold, N, f)
    ϵ_v = contagion_validity(G, E, E_threshold, D, D_threshold, N, f)
    ϵ_c = contagion_consistency(E, E_threshold, D, D_threshold, R, R_threshold, N, f)
    ϵ = max(ϵ_t, ϵ_v, ϵ_c)
    return ϵ
end
