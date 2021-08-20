function forward_nolog!(α, c, obs_density::AbstractMatrix, hmm::HMM)
    T, S = size(obs_logdensity, 1), nb_states(hmm)
    π0, P = initial_distribution(hmm), transition_matrix(hmm)
    for i = 1:S
        α[1, i] = π0[i] * obs_density[1, i]
    end
    c[1] = sum(@view α[1, :])  # scaling
    for i = 1:S
        α[1, i] /= c[1]
    end
    for t = 1:T-1
        for j = 1:S
            α[t+1, j] = sum(α[t, i] * P[i, j]) * obs_density[t+1, j]
        end
        c[t+1] = sum(@view α[t+1, :])  # scaling
        for j = 1:S
            α[t+1, j] /= c[t+1]
        end
    end
    for t = 1:T
        if all_zeros(@view α[t, :])
            throw(OverflowError("Probabilities are too small in forward step."))
        end
    end
    return nothing
end

function forward_log!(logα, obs_logdensity::AbstractMatrix, hmm::HMM)
    T, S = size(obs_logdensity, 1), nb_states(hmm)
    logπ0, logP = log.(initial_distribution(hmm)), log.(transition_matrix(hmm))
    for i = 1:S
        logα[1, i] = logπ0[i] + obs_logdensity[1, i]
    end
    for t = 1:T-1
        for j = 1:S
            logα[t+1, j] = logsumexp(logα[t, :] + logP[:, j]) + obs_logdensity[t+1, j]  # TODO
        end
    end
    for t = 1:T
        if all_nan(@view logα[t, :])
            throw(OverflowError("Log probabilities are too small in forward step."))
        end
    end
    return nothing
end

function backward_nolog!(β, c, obs_density::AbstractMatrix, hmm::HMM)
    T, S = size(obs_logdensity, 1), nb_states(hmm)
    P = transition_matrix(hmm)

    for i = 1:S
        β[T, i] = 0.0
    end
    for t = T-1:-1:1
        for i = 1:S
            β[t, i] = sum(P[i, j] * obs_density[t+1, j] * β[t+1, j] for j = 1:S) / c[t]
        end
    end
    for t = 1:T
        if all_zeros(@view β[t, :])
            throw(OverflowError("Log probabilities are too small in backward step."))
        end
    end
    return nothing
end


function backward_log!(logβ, obs_logdensity::AbstractMatrix, hmm::HMM)
    T, S = size(obs_logdensity, 1), nb_states(hmm)
    logP = log.(transition_matrix(hmm))

    for i = 1:S
        logβ[T, i] = 0.0
    end
    for t = T-1:-1:1
        for i = 1:S
            logβ[t, i] = logsumexp(logP[i, :] .+ obs_logdensity[t+1, :] .+ logβ[t+1, :])  # TODO
        end
    end
    for t = 1:T
        if all_nan(@view logβ[t, :])
            throw(OverflowError("Log probabilities are too small in backward step."))
        end
    end
    return nothing
end

function forward_backward_nolog!(α, β, c, γ, ξ, obs_density::AbstractMatrix, hmm::HMM)
    T, S = size(obs_logdensity, 1), nb_states(hmm)

    forward_nolog!(α, c, hmm, obs_density)
    backward_nolog!(β, c, hmm, obs_density)

    for t = 1:T
        for i = 1:S
            γ[t, i] = α[t, i] * β[t, i]
        end
        sumγ = sum(@view γ[t, :])
        for i = 1:S
            γ[t, i] /= sumγ
        end
    end

    for t = 1:T-1
        for i = 1:S, j = 1:S
            ξ[t, i, j] = α[t, i] * P[i, j] * obs_density[t+1, j] * β[t+1, j]
        end
        sumξ = sum(@view ξ[t, :, :])
        for i = 1:S, j = 1:S
            ξ[t, i, j] /= sumξ
        end
    end

    logL = sum(log.(c))  # TODO: @view

    return log(L)
end


function forward_backward_log!(logα, logβ, logγ, logξ, obs_logdensity::AbstractMatrix, hmm::HMM)
    T, S = size(obs_logdensity, 1), nb_states(hmm.transitions)
    logP = log.(transition_matrix(hmm))

    forward_log!(logα, obs_logdensity, hmm)
    backward_log!(logβ, obs_logdensity, hmm)

    for t = 1:T
        for i = 1:S
            logγ[t, i] = logα[t, i] + logβ[t, i]
        end
        logsumγ = logsumexp(@view logγ[t, :])
        for i = 1:S
            logγ[t, i] -= logsumγ
        end
    end

    for t = 1:T-1
        for i = 1:S, j = 1:S
            logξ[t, i, j] = logα[t, i] + logP[i, j] + obs_logdensity[t+1, j] + logβ[t+1, j]
        end
        logsumξ = logsumexp(@view logξ[t, :, :])
        for i = 1:S, j = 1:S
            logξ[t, i, j] -= logsumξ
        end
    end

    logL = logsumexp(logα[T, :])

    return logL
end

## Likelihood of observations

function update_obs_density!(obs_density::AbstractMatrix, hmm::HMM, observations::AbstractVector)
    T, S = length(observations), nb_states(hmm)
    for t = 1:T
        for s = 1:S
            obs_density[t, s] = density(emission(hmm, s), observations[t])
        end
    end
    for t = 1:T
        if all_zeros(@view obs_logdensity[t, :])
            throw(OverflowError("Densities are too small for observations."))
        end
    end
end

function update_obs_logdensity!(
    obs_logdensity::AbstractMatrix,
    hmm::HMM,
    observations::AbstractVector,
)
    T, S = length(observations), nb_states(hmm)
    for t = 1:T
        for s = 1:S
            obs_logdensity[t, s] = logdensity(emission(hmm, s), observations[t])
        end
    end
    for t = 1:T
        if all_minus_inf(@view obs_logdensity[t, :])
            throw(OverflowError("Log-densities are too small for observations."))
        end
    end
end

function baum_welch_nolog!(
    α,
    β,
    c,
    γ,
    ξ,
    obs_density::AbstractMatrix,
    hmm::HMM{Tr,Em},
    observations::AbstractVector;
    iterations,
) where {Tr,Em}
    logL_evolution = Float64[]
    for _ = 1:iterations
        update_obs_density!(obs_density, hmm, observations)
        logL = forward_backward_nolog!(α, β, c, γ, ξ, obs_density, hmm)
        push!(logL_evolution, logL)
        new_transitions = fit(Tr, γ = γ, ξ = ξ)
        new_emissions = [fit(Em, observations, γ[:, s]) for s = 1:nb_states(hmm)]  # TODO: @view
        hmm = HMM(new_transitions, new_emissions)
    end
    return hmm, logL_evolution
end

function baum_welch_log!(
    logα,
    logβ,
    logγ,
    logξ,
    obs_logdensity::AbstractMatrix,
    hmm::HMM{Tr,Em},
    observations::AbstractVector;
    iterations,
) where {Tr,Em}
    logL_evolution = Float64[]
    for _ = 1:iterations
        update_obs_logdensity!(obs_logdensity, hmm, observations)
        logL = forward_backward_log!(logα, logβ, logγ, logξ, obs_logdensity, hmm)
        push!(logL_evolution, logL)
        new_transitions = fit(Tr, exp.(logγ), exp.(logξ))
        new_emissions = [fit(Em, observations, exp.(logγ[:, s])) for s = 1:nb_states(hmm)]  # TODO: @view
        hmm = HMM(new_transitions, new_emissions)
    end
    return hmm, logL_evolution
end

## Non mutating version

function baum_welch(hmm::HMM, observations::AbstractVector; iterations, log = true)
    T, S = length(observations), nb_states(hmm)
    if log
        logα = Matrix{Float64}(undef, T, S)
        logβ = Matrix{Float64}(undef, T, S)
        logγ = Matrix{Float64}(undef, T, S)
        logξ = Array{Float64,3}(undef, T - 1, S, S)
        obs_logdensity = Matrix{Float64}(undef, T, S)
        hmm_est, logL_evolution = baum_welch_log!(
            logα,
            logβ,
            logγ,
            logξ,
            obs_logdensity,
            hmm,
            observations;
            iterations = iterations,
        )
    else
        α = Matrix{Float64}(undef, T, S)
        β = Matrix{Float64}(undef, T, S)
        c = Vector{Float64}(undef, T, S)
        γ = Matrix{Float64}(undef, T, S)
        ξ = Array{Float64,3}(undef, T - 1, S, S)
        obs_density = Matrix{Float64}(undef, T, S)
        hmm_est, logL_evolution = baum_welch_nolog!(
            α,
            β,
            c,
            γ,
            ξ,
            obs_density,
            hmm,
            observations;
            iterations = iterations,
        )
    end
    return hmm_est, logL_evolution
end
