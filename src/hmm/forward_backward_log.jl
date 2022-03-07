function forward_log!(logα::AbstractMatrix, obs_logdensity::AbstractMatrix, hmm::HMM)
    T, S = size(obs_logdensity, 1), nb_states(hmm)
    logπ0, logP = log.(initial_distribution(hmm)), log.(transition_matrix(hmm))
    for i in 1:S
        logα[1, i] = logπ0[i] + obs_logdensity[1, i]
    end
    for t in 1:(T - 1)
        for j in 1:S
            logα[t + 1, j] = logsumexp(logα[t, :] + logP[:, j]) + obs_logdensity[t + 1, j]  # TODO
        end
    end
    for t in 1:T
        if all_nan(@view logα[t, :])
            throw(OverflowError("Log probabilities are too small in forward step."))
        end
    end
    return nothing
end

function backward_log!(logβ::AbstractMatrix, obs_logdensity::AbstractMatrix, hmm::HMM)
    T, S = size(obs_logdensity, 1), nb_states(hmm)
    logP = log.(transition_matrix(hmm))

    for i in 1:S
        logβ[T, i] = 0.0
    end
    for t in (T - 1):-1:1
        for i in 1:S
            logβ[t, i] = logsumexp(logP[i, :] .+ obs_logdensity[t + 1, :] .+ logβ[t + 1, :])  # TODO
        end
    end
    for t in 1:T
        if all_nan(@view logβ[t, :])
            throw(OverflowError("Log probabilities are too small in backward step."))
        end
    end
    return nothing
end

function forward_backward_log!(
    logα::AbstractMatrix,
    logβ::AbstractMatrix,
    logγ::AbstractMatrix,
    logξ::AbstractArray{<:Real,3},
    obs_logdensity::AbstractMatrix,
    hmm::HMM,
)
    T, S = size(obs_logdensity, 1), nb_states(hmm.transitions)
    logP = log.(transition_matrix(hmm))

    forward_log!(logα, obs_logdensity, hmm)
    backward_log!(logβ, obs_logdensity, hmm)

    for t in 1:T
        for i in 1:S
            logγ[t, i] = logα[t, i] + logβ[t, i]
        end
        logsumγ = logsumexp(@view logγ[t, :])
        for i in 1:S
            logγ[t, i] -= logsumγ
        end
    end

    for t in 1:(T - 1)
        for i in 1:S, j in 1:S
            logξ[t, i, j] =
                logα[t, i] + logP[i, j] + obs_logdensity[t + 1, j] + logβ[t + 1, j]
        end
        logsumξ = logsumexp(@view logξ[t, :, :])
        for i in 1:S, j in 1:S
            logξ[t, i, j] -= logsumξ
        end
    end

    logL = logsumexp(logα[T, :])

    return logL
end

## Likelihood of observations

function update_obs_logdensity!(
    obs_logdensity::AbstractMatrix, observations::AbstractVector, hmm::HMM
)
    T, S = length(observations), nb_states(hmm)
    for t in 1:T
        for s in 1:S
            obs_logdensity[t, s] = logdensityof(emission(hmm, s), observations[t])
        end
    end
    for t in 1:T
        if all_minus_inf(@view obs_logdensity[t, :])
            throw(OverflowError("Log-densities are too small for observations."))
        end
    end
end
