include("sieve.jl")

states = [[]]

function contagion_validity(G, E, E_threshold, D, D_threshold, N::Int64=1024, f::Float64=0.1)
    ϵ_o = byzantine_echo_threshold(D, D_threshold, f)
    ϵ_pcb = sieve_total_validity(G, E, E_threshold, N, f)
    ϵ_v = ϵ_pcb + (1 - ϵ_pcb)*ϵ_o
    return ϵ_v
end

function compute_infected_round_r_only(R, R_threshold, Nbar_i, Ubar_i, l, N)
    # P[W_r_i | not(W_r_i-1), s_h]
    res = 0.0
    tmp = []
    for vim1 in 1:R_threshold
        V_i_m_1 = vim1 - 1
        push!(tmp, binomial_k(V_i_m_1, R, l*(Nbar_i-Ubar_i)/N))
    end
    denom = sum(tmp)
    for vim1 in 1:R_threshold
        V_i_m_1 = vim1-1
        prob1 = l*(Nbar_i-Ubar_i)/N
        p1 = binomial_k(V_i_m_1, R, l*(Nbar_i-Ubar_i)/N)/denom
        p2 = 0.0
        for vi in R_threshold+1:R+1
            V_i = vi-1
            p2 = p2 + binomial_k(V_i-V_i_m_1, R-V_i_m_1, (l*Ubar_i)/(N-l*(Nbar_i-Ubar_i)))
        end
        res = res + (p1*p2)
    end
    return res
end

function compute_nextStep(R, R_threshold, l, N, Nbar_r_i, Ubar_i, Nbar_next, Ubar_next)::Float64
    # Compute probability of going to the next state (N_i+1,U_i+1) given (N_i,U_i).
    # P[N_i+1, U_i+1 | N_i, U_i]
    pw = compute_infected_round_r_only(R, R_threshold, Nbar_r_i, Ubar_i, l, N)
    if pw > 1
        pw = 1
    end
    next = binomial_k(Ubar_next, N-Nbar_r_i, pw)
    return next
end

function compute_gamma(N, R, l, K, S, R_threshold, γ_end, N_r_i, old_U)::Float64
    # Compute the probability the Threshold Contagion ends with exactly γ_end
    # after K rounds. P[γ(N,R,l,K,S,R_threshold)=γ_end]

    # If initial infection S is bigger than threshold, probability is 0.
    if S > γ_end
        return 0.0
    end
    # If initial infection S is equal to threshold, probability equals the
    # the probability of going in state (N^r_i=γ_end, U^r_i=0).
    if S == γ_end
        # Pr[Ni+1=S, Ui+1=0 | Ni=S, Ui=S]
        r = compute_nextStep(R, R_threshold, l, N, S, S, γ_end, 0)
        return r
    end
    # Initialize 3 dimensional matrix of size γ_end, γ_end-S, K.
    # Value [x,y,z] contains probability of reaching state (γ_end, 0) from state (x-1,y-1) at round z.
    a = zeros(γ_end+1, γ_end-S+1, K)
    # State (γ_end, 0) at last round is 1
    a[γ_end+1, 1, 1] = 1.0
    res = 0.0
    for k in 1:K
        for n in γ_end+1:-1:((K-k+1)*S)+1
            act_N_r_i = n - 1
            # Case where U^r_i = 0
            if n < γ_end+1
                # If last round, absorbing state
                if k == 1
                    a[n, 1, k] = 0.0
                else
                    # Condition of equation 32
                    if n+S > N
                        a[n, 1, k] = a[n, 1, k-1]
                    # If infecting S new nodes at beginning of the round
                    # overflows the threshold, absorbing state
                    elseif n+S > γ_end+1
                        a[n, 1, k] = 0.0
                    else
                        a[n, 1, k] = a[n+S, S, k-1]
                    end
                end
            end
            # If we are at state (N^r_i=S, U^r_i) at round K (starting round),
            # compute probability of ending in state (N^r_i=γ_end, U^r_i=0)
            if act_N_r_i == S && k == K
                tmp_res = 0.0
                # Each possible incrementation of infected nodes (0 to γ-Ni)
                for add in 1:(γ_end-act_N_r_i)+1
                    new_infected = add-1
                    tmp = a[n+add-1, add, k]
                    if tmp > 0
                        # Pr[Ni+1=n+add, Ui+1=add | Ni=n, Ui=u]
                        p = compute_nextStep(R, R_threshold, l, N, S, S, act_N_r_i+new_infected, new_infected)
                        tmp_res = tmp_res + p*tmp
                    end
                end
                res = tmp_res
            else
                for u in 2:n-S
                    act_U_r_i = u-1
                    tmp_res = 0.0
                    # Each possible incrementation of infected nodes (0 to γ_end-Ni)
                    for add in 1:(γ_end-act_N_r_i-(K-k)*S)+1
                        new_infected = add-1
                        tmp = a[n+add-1, add, k]
                        if tmp > 0
                            # Pr[Ni+1=n+add, Ui+1=add | Ni=n, Ui=u]
                            p = compute_nextStep(R, R_threshold,l, N, act_N_r_i, act_U_r_i, act_N_r_i+new_infected, new_infected)
                            tmp_res = tmp_res + p*tmp
                        end
                    end
                    a[n, u, k] = tmp_res
                end
            end
        end
    end
    return res
end

function specific_enough_ready(D, D_threshold, γ_barPlus, N)
    # Compute probability that one specific process eventually collects enough Ready(m)
    # to deliver m.
    μ_tilde = 0.0
    for D_bar in D_threshold:D
        μ_tilde = μ_tilde + binomial_k(D_bar, D, γ_barPlus/N)
    end
    if μ_tilde > 1
        return 1.0
    end
    return μ_tilde
end

function any_enough_ready(N, C, D, D_threshold, R, R_threshold)
    # Compute probability any correct process eventually collects enough Ready(m)
    # to deliver m.
    μ = 0.0
    counter = 5
    for γ_barPlus in N-C:N
        p_g = compute_gamma(N, R, 1, 1, N-C, R_threshold, γ_barPlus, 0, 0)
        μ = μ + (1-(1-specific_enough_ready(D, D_threshold, γ_barPlus, N))^C)*p_g
        if p_g == 0
            counter = counter - 1
        end
        if counter == 0
            break
        end
    end
    return μ
end

function contagion_consistency(E, E_threshold, D, D_threshold, R, R_threshold, N::Int64=1024, f::Float64=0.1)
    C = floor(Int, (1-f)*N)
    μ = any_enough_ready(N, C, D, D_threshold, R, R_threshold)
    ϵ_pcb = sieve_consistency(E, E_threshold, N, f)
    ϵ_c = ϵ_pcb + (1 - ϵ_pcb)*μ
    return ϵ_c
end

function alpha_minus(γ, N, D, D_threshold)
    res = 0.0
    for D_bar in D_threshold:D
        res = res + binomial_k(D_bar, D, γ/N)
    end
    return res
end

function alpha_plus(γ, N, C, D, D_threshold)
    res = 0.0
    for D_bar in D_threshold:D
        res = res + binomial_k(D_bar, D, (γ+(N-C))/N)
    end
    return res
end

function compute_alpha(γ, N, C, D, D_threshold)
    # Upper bound of probability of totality being compromised given γ correct
    # processes eventually ready for m.
    α_minus = alpha_minus(γ, N, D, D_threshold)
    α_plus = alpha_plus(γ, N, C, D, D_threshold)
    res = 1 - (α_minus^C) - (1 - α_plus)^C
    return res
end

function epsilon_b(N, f, C, D, D_threshold, R, R_threshold)
    # Probability totality being compromised. The adversary plays up to C rounds,
    # and up to C correct processes can be eventually ready for m.
    ϵ_b = 0.0
    for n in 1:C
        for γ_n in 1:C
            # If n is too small, the ready loop won't start
            if n < R_threshold
                # End of round n there are n infected nodes
                if γ_n == n
                    p_γ = compute_gamma(C, R, 1-f, n, 1, R_threshold, γ_n, 0, 0)
                    α = compute_alpha(γ_n, N, C, D, D_threshold)
                    ϵ_b = ϵ_b + (p_γ*α)
                else
                    p_γ = 0.0
                end
            else
                # At round n there are at least n infected nodes (1 per round beginning)
                if γ_n >= n
                    p_γ = compute_gamma(C, R, 1-f, n, 1, R_threshold, γ_n, 0, 0)
                    α = compute_alpha(γ_n, N, C, D, D_threshold)
                    ϵ_b = ϵ_b + (p_γ*α)
                else
                    p_γ = 0.0
                end
            end
        end
    end
    return ϵ_b
end

function contagion_totality(E, E_threshold, D, D_threshold, R, R_threshold, N::Int64=1024, f::Float64=0.1)
    C = floor(Int, (1-f)*N)
    ϵ_pcb = sieve_consistency(E, E_threshold, N, f)
    μ = any_enough_ready(N, C, D, D_threshold, R, R_threshold)
    ϵ_b = epsilon_b(N, f, C, D, D_threshold, R, R_threshold)
    ϵ_t = ϵ_pcb + μ + ϵ_b
end

function runPwi()
    N = 8192;
    R = 128;
    R_threshold = 48;
    Ni = 2000
    f = 0.15;

    for u in 0:Ni
        v = compute_infected_round_r_only(R, R_threshold, Ni, u, 1, N)
        println("($Ni/$u): $v")
    end
end
