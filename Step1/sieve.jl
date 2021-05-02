include("murmur.jl")

function binomial_k(k, n, p)
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
    res = log_binomial_coefficient(n, k) + log(p)*k + log(1-p)*(n-k)
    return Float64(exp(res))
end

function byzantine_echo_threshold(E, E_threshold, f)
    # Compute probability of having more than E-E_threshold Byzantine
    # processes in the sampling set E.
    ϵ_o = 0.0
    for F_bar in E-E_threshold+1:E
        ϵ_o = ϵ_o + binomial_k(F_bar, E, f)
    end
    return ϵ_o
end

function sieve_total_validity(G, E, E_threshold, N::Int64=1024, f::Float64=0.1)
    C = floor(Int, (1-f)*N)
    ϵ_t = murmur_totality(G, N, f)
    ϵ_o = byzantine_echo_threshold(E, E_threshold, f)
    ϵ_v = ϵ_t + (1-ϵ_t)*(1-(1-ϵ_o)^C)
    return ϵ_v
end

function alpha_function(M, F_barPi, C, E, E_threshold)
    # α(M, F_barPi) as described at page 95.
    if (E-F_barPi) == 0 || (M/C) == 0
        α = 0
    else
        α = ((exp(1)*(E-F_barPi)*(M/C))/(E_threshold-F_barPi))^(E_threshold-F_barPi)
    end
    return α
end

function beta_function(M, F_barPi, C, E)
    # β(M, F_barPi) as described at page 95.
    β = exp(-(E-F_barPi)*(M/C))
    return β
end

function simple_psi_two_params(M, F_barPi, C, E, E_threshold)
    # Compute probability that at least L processes pb.Deliver a message m,
    # before π delivers m, given F_barPi Byzantine processes in π's first echo sample.
    if E_threshold < F_barPi
        psi_condition = 0
    else
        psi_condition = (((E_threshold-F_barPi)-sqrt(E_threshold-F_barPi))/(E-F_barPi))
    end
    if (M/C) <= psi_condition
        ψ = alpha_function(M, F_barPi, C, E, E_threshold)*beta_function(M, F_barPi, C, E)
    else
        ψ = 1
    end
    return ψ
end

function simple_psi(L, C, E, E_threshold, f)
    # Compute probability that at least L processes pb.Deliver a message m,
    # before π delivers m.
    ψ = 0.0
    for F_barPi in 0:E_threshold
        ψ = ψ + (1 - (1-simple_psi_two_params(L, F_barPi, C, E, E_threshold))^floor(Int, C/L) * (1-simple_psi_two_params(mod(C, L), F_barPi, C, E, E_threshold)))*binomial_k(F_barPi, E, f)
    end
    return ψ
end

function psi_tilde(L, C, E, E_threshold, f)
    # Compute probability that exactly L processes pb.Deliver a message m,
    # before π delivers m.
    if L == C
        ψ_tilde = 1
    elseif L == 1
        ψ_tilde = 1-(1-simple_psi(L, C, E, E_threshold, f))^C
    else
        ψ_tilde = (1-(1-simple_psi(L, C, E, E_threshold, f))^C) - (1-(1-simple_psi(L-1, C, E, E_threshold, f))^C)
    end
    return ψ_tilde
end

function deliverM_after_kProcesses_given_byzantinePopulation(E, E_threshold, F_barPi, k, C)
    # P[A^k_m[π] | F_π]
    res = 0.0
    for echoes in E_threshold-F_barPi:E-F_barPi
        res = res + binomial_k(echoes, E-F_barPi, k/C)
    end
    # Due to rounding, sum can be slightly bigger than 1
    if res > 1
        res = 1.0
    end
    return Float64(res)
end

function phi_plus(C, f, E, E_threshold, N_barH)
    # ∑ ϕ(N_H) * P[F_π^+ | N_H]
    sub_res = []
    # Consider non poisoned systems
    for F_barPi in 0:E_threshold-1
        P_F_barPi = binomial_k(F_barPi, E, f)
        tmp = (deliverM_after_kProcesses_given_byzantinePopulation(E, E_threshold, F_barPi, N_barH, C)-deliverM_after_kProcesses_given_byzantinePopulation(E, E_threshold, F_barPi, N_barH-1, C))*P_F_barPi
        # Due to rounding the difference can be a very small negative value.
        if tmp < 0
            tmp = 0.0
        end
        push!(sub_res, tmp)
    end
    denom = sum(sub_res)
    res = 0.0
    # Consider non poisoned systems
    if denom == 0
        res = 0
    else
        for F_barPi in 0:E_threshold-1
            res = res + simple_phi(F_barPi, N_barH, C, E, E_threshold)*(sub_res[F_barPi+1]/denom)
        end
    end
    # Due to rounding, sum can be slightly bigger than 1
    if res > 1
        res = 1.0
    end
    return res
end

function phi_minus(C, f, E, E_threshold, N_barH)
    # ∑ ϕ(N_H) * P[F_π^- | N_H]
    sub_res = []
    # Consider non poisoned systems
    for F_barPi in 0:E_threshold-1
        P_F_barPi = binomial_k(F_barPi, E, f)
        tmp = deliverM_after_kProcesses_given_byzantinePopulation(E, E_threshold, F_barPi, N_barH-1, C)*P_F_barPi
        push!(sub_res, tmp)
    end
    denom = 1-sum(sub_res)
    res = 0.0
    # Consider non poisoned systems
    for F_barPi in 0:E_threshold-1
        P_F_barPi = binomial_k(F_barPi, E, f)
        res = res + simple_phi(F_barPi, N_barH, C, E, E_threshold) * (((1-sub_res[F_barPi+1])*P_F_barPi)/denom)
    end
    return res
end

function simple_phi(F_barPi, N_barHbar, C, E, E_threshold)
    # Compute bound of the probability that the a correct process
    # delivers any message m ≂̸ H, given N_barH the number of processes
    # that pb.Delivered H and F_barPi the number of Byzantine processes
    # in the first echo sample.
    if E_threshold < F_barPi
        phi_condition = 0
    else
        phi_condition = ((E_threshold-F_barPi)-sqrt(E_threshold-F_barPi))/(E-F_barPi)
    end
    if (C-N_barHbar)/C <= phi_condition
        ϕ = alpha_function(C-N_barHbar, F_barPi, C, E, E_threshold)*beta_function(C-N_barHbar, F_barPi, C, E)
    else
        ϕ = 1
    end
    return ϕ
end

function phi_tilde(N_barH, C, E, E_threshold, f)
    # Compute bound of the probability that the adversary compromises
    # consistency, given N_barH processes have pb.Delivered H.
    ϕ_plus = phi_plus(C, f, E, E_threshold, N_barH)
    ϕ_minus = phi_minus(C, f, E, E_threshold, N_barH)
    # Due to rounding, the value can be a very small negative number
    if ϕ_minus < 0
        ϕ_minus = 0.0
    end
    # Due to rounding, sum can be slightly bigger than 1
    if ϕ_minus > 1
        ϕ_minus = 1.0
    end
    ϕ = (1 - (1 - ϕ_plus)*(1-ϕ_minus)^(C-1))
    return ϕ
end

function epsilon_p(C, f, E, E_threshold)
    # Compute probability of system being poisoned.
    tmp = 0.0
    for F_bar in E_threshold:E
        tmp = tmp + binomial_k(F_bar, E, f)
    end
    ϵ_p = 1 - (1 - tmp)^C
    return ϵ_p
end

function sieve_consistency(E, E_threshold, N::Int64=1024, f::Float64=0.1)
    C = floor(Int, (1-f)*N)
    ϵ_c = epsilon_p(C, f, E, E_threshold)
    for L in 1:C
        ϵ_c = ϵ_c + phi_tilde(L, C, E, E_threshold, f)*psi_tilde(L, C, E, E_threshold, f)
    end
    return ϵ_c
end
