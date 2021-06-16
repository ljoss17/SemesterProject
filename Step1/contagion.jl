include("sieve.jl")

mutable struct InfectionProbabilities
    N::Int64
    l::Float64
    R::Int64
    R_thr::Int64
    prob::Array{Float128,2}
end

pwi = InfectionProbabilities(0, 0.0, 0, 0, zeros(1,1))

function contagion_validity(G, E, E_threshold, D, D_threshold, N::Int64=1024, f::Float64=0.1)::Float128
    ϵ_o::Float128 = sum_binomial(D-D_threshold+1, D, D, f)
    ϵ_pcb::Float128 = sieve_total_validity(G, E, E_threshold, N, f)
    ϵ_v::Float128 = ϵ_pcb + (1 - ϵ_pcb)*ϵ_o
    return check_rounding(ϵ_v)
end

function compute_infected_round_r_only(R, R_threshold, Nbar_i, Ubar_i, l, N)::Float128
    # P[W_r_i+1 | not(W_r_i-1), s_h]
    denom::Float128 = sum_binomial(0, R_threshold-1, R, l*(Nbar_i-Ubar_i)/N)
    vals::Array{Float128, 1} = zeros(R_threshold)
    for vim1 in 1:R_threshold
        Vim1::Int64 = vim1-1
        p_Wim1::Float128 = l*(Nbar_i-Ubar_i)/N
        p_Vim1_nWi::Float128 = binomial_k(R, p_Wim1, Vim1)/denom
        tmp_vals::Array{Float128, 1} = zeros(R-R_threshold+1)
        for vi in R_threshold+1:R+1
            Vi::Int64 = vi-1
            tmp_vals[vi-R_threshold] = binomial_k(R-Vim1, (l*Ubar_i)/(N-l*(Nbar_i-Ubar_i)), Vi-Vim1)
        end
        p_Vi_Vim1_nWi::Float128 = kahan_summation(tmp_vals)
        vals[vim1] = (p_Vim1_nWi*p_Vi_Vim1_nWi)
    end
    res::Float128 = kahan_summation(vals)
    return res
end

function compute_nextStep(R, R_threshold, l, N, Nbar_r_i, Ubar_i, Nbar_next, Ubar_next)::Float128
    # Compute probability of going to the next state (N_i+1,U_i+1) given (N_i,U_i).
    # P[N_i+1, U_i+1 | N_i, U_i]
    pw::Float128 = pwi.prob[Nbar_r_i+1, Ubar_i+1]
    next::Float128 = binomial_k(N-Nbar_r_i, pw, Ubar_next)
    return next
end

function compute_all_steps(N, l, R, R_threshold, ni_min, ni_max)
    pwi.prob = zeros(ni_max+1, ni_max+1)
    # If there are less than R_threshold infected nodes, the probability of a node
    # being infected is 0.
    for ni in R_threshold+1:ni_max+1
        Ni::Int64 = ni-1
        for ui in 1:ni
            Ui::Int64 = ui-1
            Wip1_nWi::Float128 = compute_infected_round_r_only(R, R_threshold, Ni, Ui, l, N)
            # Due to rounding, sum can be slightly bigger than 1
            pwi.prob[ni, ui] = min(1.0, Wip1_nWi)
        end
    end
    global pwi.N = N
    global pwi.l = l
    global pwi.R = R
    global pwi.R_thr = R_threshold
end

function threshold_contagion_single_round(N, R, l, K, S, R_threshold, γ_end)
    if (pwi.N != N) || (pwi.l != l) || (pwi.R != R) || (pwi.R_thr != R_threshold)
        compute_all_steps(N, l, R, R_threshold, 0, N)
    end
    dist = zeros(γ_end+1, γ_end+1)
    # Probability of going from (S,S) to (S, 0)
    for u in 0:γ_end-S
        n::Int64 = S+u
        dist[n+1, u+1] = compute_nextStep(R, R_threshold, l, N, S, S, n, u)
    end
    for N_i in S+1:γ_end
        n_i::Int64 = N_i+1
        for U_i in N_i-S:-1:0
            u_i::Int64 = U_i+1
            previous_n::Int64 = N_i-U_i
            for i in 1:previous_n
                t1 = dist[previous_n+1, i+1]
                if t1 != 0
                    t2::Float128 = compute_nextStep(R, R_threshold, l, N, previous_n, i, N_i, U_i)
                    dist[n_i, u_i] = dist[n_i, u_i] + t1*t2
                end
            end
        end
    end
    return dist[:,1]
end

function threshold_contagion_multiple_rounds(N, R, l, K, S, R_threshold, γ_end, dist_old)
    if (pwi.N != N) || (pwi.l != l) || (pwi.R != R) || (pwi.R_thr != R_threshold)
        compute_all_steps(N, l, R, R_threshold, 0, N)
    end
    dist = copy(dist_old)
    if K == 1
        for u in 0:γ_end-S
            n::Int64 = S+u
            dist[n+1, u+1] = compute_nextStep(R, R_threshold, l, N, S, S, n, u)
        end
    end
    for N_i in K:γ_end
        n_i::Int64 = N_i+1
        for U_i in N_i-S:-1:0
            u_i::Int64 = U_i+1
            previous_n::Int64 = N_i-U_i
            for i in 1:previous_n
                t1 = dist[previous_n+1, i+1]
                if t1 != 0
                    t2::Float128 = compute_nextStep(R, R_threshold, l, N, previous_n, i, N_i, U_i)
                    dist[n_i, u_i] = dist[n_i, u_i] + t1*t2
                end
            end
        end
    end
    return dist[:,1]
end

function contagion_threshold(N, R, l, K, S, R_threshold, γ_end)
    dist_old = zeros(γ_end+1, γ_end+1)
    all_dists = fill(Float128[], 1, K)
    for k in 1:K
        dist = threshold_contagion_multiple_rounds(N, R, l, k, S, R_threshold, γ_end, dist_old)
        all_dists[k] = copy(dist)
        dist_old = zeros(γ_end+1, γ_end+1)
        dist_old[3:γ_end+1, 2] = dist[2:γ_end]
        # If infecting S new nodes is results in more than N, then keep N infected.
        # So the probability of being in state (γ_end, 0) at round k+1 is equal to
        # the probability of being in state (γ_end, 0) at round k, plus the contagion threshold computations.
        dist_old[γ_end+1, 1] = dist[γ_end+1]
    end
    return all_dists
end

function compute_gamma(N, R, l, K, S, R_threshold, γ_end)::Float64
    # This function is not used.

    # Compute the probability the Threshold Contagion ends with exactly γ_end
    # after K rounds. P[γ(N,R,l,K,S,R_threshold)=γ_end]

    # If initial infection S, times the number of rounds, is bigger than threshold,
    # probability is 0 (more than γ nodes will be infected).
    if (S*K) > γ_end
        return 0.0
    end

    # If at the end of all rounds there are less than R_threshold infected nodes
    # and gamma is bigger than this number, probability is 0 (only S*K nodes will be infected).
    if ((S*K) < R_threshold) && (γ_end > (S*K))
        return 0.0
    end

    # If at the end of all rounds there are less than R_threshold infected nodes
    # exactly S*K nodes will be infected. So if γ == S*K, probability is 1.
    if ((S*K) < R_threshold) && (γ_end == (S*K))
        return 1.0
    end

    # If the probabilities for a node to be infected are not for the current
    # system, compute them before starting the Contagion Threshold.
    if (pwi.N != N) || (pwi.l != l) || (pwi.R != R) || (pwi.R_thr != R_threshold)
        compute_all_steps(N, l, R, R_threshold, 0, N)
    end

    # If initial infection S is equal to threshold, probability equals the
    # the probability of going in state (N^r_i=γ_end, U^r_i=0).
    # This only takes into account if K = 1, if K > 1 and S == γ, first if returns 0.
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

function specific_enough_ready(D, D_threshold, γ_barPlus, N)::Float128
    # Compute probability that one specific process eventually collects enough Ready(m)
    # to deliver m.
    μ_tilde::Float128 = sum_binomial(D_threshold, D, D, γ_barPlus/N)
    return check_rounding(μ_tilde)
end

function any_enough_ready(N, C, D, D_threshold, R, R_threshold)::Float128
    # Compute probability any correct process eventually collects enough Ready(m)
    # to deliver m.
    dists = contagion_threshold(N, R, 1, 1, N-C, R_threshold, N)
    γs = dists[1]
    vals = zeros(C+1)
    for γ_barPlus in N-C:N
        #p_g = compute_gamma(N, R, 1, 1, N-C, R_threshold, γ_barPlus)
        p_g::Float128 = γs[γ_barPlus+1]
        #μ = μ + (1-(1-specific_enough_ready(D, D_threshold, γ_barPlus, N))^C)*p_g
        vals[γ_barPlus-N+C+1] = (1-(1-specific_enough_ready(D, D_threshold, γ_barPlus, N))^C)*p_g
    end
    μ::Float128 = kahan_summation(vals)
    return μ
end

function contagion_consistency(E, E_threshold, D, D_threshold, R, R_threshold, N::Int64=1024, f::Float64=0.1)::Float128
    C = floor(Int, (1-f)*N)
    μ::Float128 = any_enough_ready(N, C, D, D_threshold, R, R_threshold)
    ϵ_pcb::Float128 = sieve_consistency(E, E_threshold, N, f)
    ϵ_c::Float128 = ϵ_pcb + (1 - ϵ_pcb)*μ
    return check_rounding(ϵ_c)
end

function alpha_minus(γ, N, D, D_threshold)::Float128
    res::Float128 = sum_binomial(D_threshold, D, D, γ/N)
    return res
end

function alpha_plus(γ, N, C, D, D_threshold)::Float128
    res::Float128 = sum_binomial(D_threshold, D, D, (γ+(N-C))/N)
    return res
end

function compute_alpha(γ, N, C, D, D_threshold)::Float128
    # Upper bound of probability of totality being compromised given γ correct
    # processes eventually ready for m.
    α_minus::Float128 = alpha_minus(γ, N, D, D_threshold)
    α_plus::Float128 = alpha_plus(γ, N, C, D, D_threshold)
    # α_minus and α_plus are lower bounds, and since α is an upper bounds,
    # the worst case is when α_minus and/or α_plus are 0.0
    P_A_γ::Float128 = max(α_minus^C, 0.0)
    P_Atilde_γ::Float128 = max((1 - α_plus)^C, 0.0)
    res::Float128 = 1 - P_A_γ - P_Atilde_γ
    return res
end

function epsilon_b(N, f, C, D, D_threshold, R, R_threshold)::Float128
    # Probability totality being compromised. The adversary plays up to C rounds,
    # and up to C correct processes can be eventually ready for m.
    dists = contagion_threshold(C, R, 1-f, C, 1, R_threshold, C)
    vals::Array{Float128, 1} = zeros(C*C)
    p_γ::Float128 = 0.0
    α::Float128 = 0.0
    for n in 1:C
        for γ_n in 1:C
            current_γ = dists[n]
            # If n is too small, the ready loop won't start
            if n < R_threshold
                # End of round n there are n infected nodes
                if γ_n == n
                    p_γ = current_γ[γ_n+1]
                    α = compute_alpha(γ_n, N, C, D, D_threshold)
                    vals[(C*(n-1))+γ_n] = (p_γ*α)
                end
            else
                # At round n there are at least n infected nodes (1 per round beginning)
                if γ_n >= n
                    p_γ = current_γ[γ_n+1]
                    α = compute_alpha(γ_n, N, C, D, D_threshold)
                    vals[(C*(n-1))+γ_n] = (p_γ*α)
                end
            end
        end
    end
    ϵ_b::Float128 = kahan_summation(vals)
    return ϵ_b
end

function contagion_totality(E, E_threshold, D, D_threshold, R, R_threshold, N::Int64=1024, f::Float64=0.1)::Float128
    C = floor(Int, (1-f)*N)
    ϵ_pcb::Float128 = sieve_consistency(E, E_threshold, N, f)
    μ::Float128 = any_enough_ready(N, C, D, D_threshold, R, R_threshold)
    ϵ_b::Float128 = epsilon_b(N, f, C, D, D_threshold, R, R_threshold)
    ϵ_t::Float128 = ϵ_pcb + μ + ϵ_b
    return check_rounding(ϵ_t)
end
