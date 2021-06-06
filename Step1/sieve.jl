include("murmur.jl")

function check_rounding(v)::Float128
    # Due to rounding, check if probability is out of the bound (0,1)
    if v > 1.0
        return 1.0
    end
    if v < 0.0
        return 0.0
    end
    return v
end

function sum_binomial(from, to, n, p)::Float128
    # Sum binomial distribution using kahan summation algorithm.
    val::Array{Float128, 1} = zeros(to-from+1)
    for i in from:to
        val[i-from+1] = binomial_k(n, p, i)
    end
    res::Float128 = kahan_summation(val)
    return res
end

function binomial_k(n, p, k)::Float128
    # Compute probability of having value k in a
    # binomial distribution Bin[n, p]

    # If probability is 1 or 0, only two possible cases
    if p == 1
        if k == n
            return 1.0
        else
            return 0.0
        end
    end
    if p == 0
        if k == 0
            return 1.0
        else
            return 0.0
        end
    end
    res::Float128 = log_binomial_coefficient(n, k) + log(p)*k + log(1-p)*(n-k)
    return exp(res)
end

function binomial_k_taylor(n, p, k)::Float128
    # Compute probability of having value k in a
    # binomial distribution Bin[n, p]

    # If probability is 1 or 0, only two possible cases

    shift_taylor(a) = a+Taylor1(typeof(a), 5)
    t = shift_taylor(0.0)
    if p == 1
        if k == n
            return 1.0
        else
            return 0.0
        end
    end
    if p == 0
        if k == 0
            return 1.0
        else
            return 0.0
        end
    end
    res::Float128 = log_binomial_coefficient(n, k) + log(p)*k + evaluate(log(1-t), p)*(n-k)
    return exp(res)
end


function byzantine_echo_threshold(E, E_threshold, f)::Float128
    # Compute probability of having more than E-E_threshold Byzantine
    # processes in the sampling set E.
    ϵ_o::Float128 = sum_binomial(E-E_threshold+1, E, E, f)
    return check_rounding(ϵ_o)
end

function sieve_total_validity(G, E, E_threshold, N::Int64=1024, f::Float64=0.1)::Float128
    C = floor(Int, (1-f)*N)
    ϵ_t::Float128 = murmur_totality(G, N, f)
    ϵ_o::Float128 = byzantine_echo_threshold(E, E_threshold, f)
    ϵ_v::Float128 = ϵ_t + (1-ϵ_t)*(1-(1-ϵ_o)^C)
    return check_rounding(ϵ_v)
end

function alpha_function(M, F_barPi, C, E, E_threshold)::Float128
    # α(M, F_barPi) as described at page 95.
    α::Float128 = 0.0
    if (E-F_barPi) == 0 || (M/C) == 0
        α = 0.0
    else
        α = ((exp(1)*(E-F_barPi)*(M/C))/(E_threshold-F_barPi))^(E_threshold-F_barPi)
    end
    return check_rounding(α)
end

function beta_function(M, F_barPi, C, E)::Float128
    # β(M, F_barPi) as described at page 95.
    β::Float128 = exp(-(E-F_barPi)*(M/C))
    return check_rounding(β)
end

function simple_psi_two_params(M, F_barPi, C, E, E_threshold)::Float128
    # Compute probability that at least L processes pb.Deliver a message m,
    # before π delivers m, given F_barPi Byzantine processes in π's first echo sample.
    psi_condition::Float128 = 0.0
    if E_threshold >= F_barPi
        psi_condition = (((E_threshold-F_barPi)-sqrt(E_threshold-F_barPi))/(E-F_barPi))
    end
    ψ::Float128 = 1.0
    if (M/C) <= psi_condition
        ψ = alpha_function(M, F_barPi, C, E, E_threshold)*beta_function(M, F_barPi, C, E)
    end
    return check_rounding(ψ)
end

function simple_psi(L, C, E, E_threshold, f)::Float128
    # Compute probability that at least L processes pb.Deliver a message m,
    # before π delivers m.
    vals::Array{Float128, 1} = zeros(E_threshold+1)
    for F_barPi in 0:E_threshold
        vals[F_barPi+1] = (1 - (1-simple_psi_two_params(L, F_barPi, C, E, E_threshold))^floor(Int, C/L) * (1-simple_psi_two_params(mod(C, L), F_barPi, C, E, E_threshold)))*binomial_k(E, f, F_barPi)
    end
    ψ::Float128 = kahan_summation(vals)
    return check_rounding(ψ)
end

function psi_tilde(L, C, E, E_threshold, f)::Float128
    # Compute probability that exactly L processes pb.Deliver a message m,
    # before π delivers m.
    ψ_tilde::Float128 = 0.0
    if L == C
        ψ_tilde = 1.0
    elseif L == 1
        ψ_tilde = 1-(1-simple_psi(L, C, E, E_threshold, f))^C
    else
        ψ_tilde = (1-(1-simple_psi(L, C, E, E_threshold, f))^C) - (1-(1-simple_psi(L-1, C, E, E_threshold, f))^C)
    end
    return check_rounding(ψ_tilde)
end

function deliverM_after_kProcesses_given_byzantinePopulation(E, E_threshold, F_barPi, k, C, f)::Float128
    # P[A^k_m[π] | F_π]
    res::Float128 = sum_binomial(E_threshold-F_barPi, E-F_barPi, E-F_barPi, k/C)
    return check_rounding(res)
end

function deliverM_after_kProcesses(E, E_threshold, k, C, f)::Tuple{Float128, Array{Float128, 1}}
    # P[A^k_m[π]]
    vals::Array{Float128, 1} = zeros(E_threshold)
    for F_barPi in 0:E_threshold-1
        pf::Float128 = binomial_k(E, f, F_barPi)
        P_ANhk_Fπ::Float128 = deliverM_after_kProcesses_given_byzantinePopulation(E, E_threshold, F_barPi, k, C, f)
        vals[F_barPi+1] = pf*P_ANhk_Fπ
    end
    res::Float128 = kahan_summation(vals)
    return check_rounding(res), vals
end

function phi_plus(C, f, E, E_threshold, N_barH)::Float128
    # ∑ ϕ(N_H) * P[F_π^+ | N_H]
    P_ANh_H::Float128, P_ANh_H_Fπ::Array{Float128, 1},  = deliverM_after_kProcesses(E, E_threshold, N_barH, C, f)
    P_ANhm1_H::Float128, P_ANhm1_H_Fπ::Array{Float128, 1} = deliverM_after_kProcesses(E, E_threshold, N_barH-1, C, f)
    denom::Float128 = P_ANh_H-P_ANhm1_H
    res_vals::Array{Float128, 1} = zeros(E_threshold)
    if denom == 0
        println("ERROR denom is 0")
        return 0.0
    end
    for F_barPi in 0:E_threshold-1
        #pf::Float128 = binomial_k(E, f, F_barPi)
        #tmp::Float128 = (deliverM_after_kProcesses_given_byzantinePopulation(E, E_threshold, F_barPi, N_barH, C, f) - deliverM_after_kProcesses_given_byzantinePopulation(E, E_threshold, F_barPi, N_barH-1, C, f))*pf
        tmp::Float128 = P_ANh_H_Fπ[F_barPi+1] - P_ANhm1_H_Fπ[F_barPi+1]
        res_vals[F_barPi+1] = simple_phi(F_barPi, N_barH, C, E, E_threshold)*tmp
    end
    res::Float128 = kahan_summation(res_vals)/denom
    return check_rounding(res)
end

function phi_minus(C, f, E, E_threshold, N_barH)::Float128
    # ∑ ϕ(N_H) * P[F_π^- | N_H]
    P_ANhm1_H::Float128, P_ANhm1_H_Fπ::Array{Float128, 1} = deliverM_after_kProcesses(E, E_threshold, N_barH-1, C, f)
    denom::Float128 = 1-P_ANhm1_H
    vals::Array{Float128, 1} = zeros(E_threshold)
    # Consider non poisoned systems
    for F_barPi in 0:E_threshold-1
        pf::Float128 = binomial_k(E, f, F_barPi)
        #tmp::Float128 = (1-deliverM_after_kProcesses_given_byzantinePopulation(E, E_threshold, F_barPi, N_barH-1, C, f))*pf
        tmp::Float128 = pf-P_ANhm1_H_Fπ[F_barPi+1]
        vals[F_barPi+1] = simple_phi(F_barPi, N_barH, C, E, E_threshold)*tmp
    end
    res::Float128 = kahan_summation(vals)/denom
    return check_rounding(res)
end

function simple_phi(F_barPi, N_barHbar, C, E, E_threshold)::Float128
    # Compute bound of the probability that the a correct process
    # delivers any message m ≂̸ H, given N_barH the number of processes
    # that pb.Delivered H and F_barPi the number of Byzantine processes
    # in the first echo sample.
    phi_condition::Float128 = 0.0
    if E_threshold < F_barPi
        phi_condition = 0.0
    else
        phi_condition = ((E_threshold-F_barPi)-sqrt(E_threshold-F_barPi))/(E-F_barPi)
    end
    ϕ::Float128 = 1.0
    if (C-N_barHbar)/C <= phi_condition
        ϕ = alpha_function(C-N_barHbar, F_barPi, C, E, E_threshold)*beta_function(C-N_barHbar, F_barPi, C, E)
    end
    return check_rounding(ϕ)
end

function phi_tilde(N_barH, C, E, E_threshold, f)::Float128
    # Compute bound of the probability that the adversary compromises
    # consistency, given N_barH processes have pb.Delivered H.
    ϕ_plus::Float128 = phi_plus(C, f, E, E_threshold, N_barH)
    ϕ_minus::Float128 = phi_minus(C, f, E, E_threshold, N_barH)
    ϕ::Float128 = (1 - (1 - ϕ_plus)*(1-ϕ_minus)^(C-1))
    return check_rounding(ϕ)
end

function epsilon_p(C, f, E, E_threshold)::Float128
    # Compute probability of system being poisoned.
    vals::Array{Float128, 1} = zeros(E-E_threshold+1)
    for F_bar in E_threshold:E
        vals[F_bar-E_threshold+1] = binomial_k(E, f, F_bar)
    end
    tmp::Float128 = kahan_summation(vals)
    ϵ_p::Float128 = 1 - (1 - tmp)^C
    return check_rounding(ϵ_p)
end

function sieve_consistency(E, E_threshold, N::Int64=1024, f::Float64=0.1)::Float128
    C = floor(Int, (1-f)*N)
    vals::Array{Float128, 1} = zeros(C+1)
    vals[1] = epsilon_p(C, f, E, E_threshold)
    for L in 1:C
        vals[L+1] = phi_tilde(L, C, E, E_threshold, f)*psi_tilde(L, C, E, E_threshold, f)
    end
    ϵ_c::Float128 = kahan_summation(vals)
    return check_rounding(ϵ_c)
end
