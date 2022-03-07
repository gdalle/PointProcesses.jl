function forward!(
    α::AbstractMatrix, c::AbstractVector, obs_density::AbstractMatrix, hmm::HMM
)
    T, S = size(obs_density, 1), nb_states(hmm)
    π0, P = initial_distribution(hmm), transition_matrix(hmm)
    for i in 1:S
        α[1, i] = π0[i] * obs_density[1, i]
    end
    c[1] = sum(@view α[1, :])  # scaling
    for i in 1:S
        α[1, i] /= c[1]
    end
    for t in 1:(T - 1)
        for j in 1:S
            α[t + 1, j] = sum(α[t, i] * P[i, j] for i in 1:S) * obs_density[t + 1, j]
        end
        c[t + 1] = sum(@view α[t + 1, :])  # scaling
        for j in 1:S
            α[t + 1, j] /= c[t + 1]
        end
    end
    for t in 1:T
        if all_zeros(@view α[t, :])
            throw(OverflowError("Probabilities are too small in forward step."))
        end
    end
    return nothing
end

function backward!(
    β::AbstractMatrix, c::AbstractVector, obs_density::AbstractMatrix, hmm::HMM
)
    T, S = size(obs_density, 1), nb_states(hmm)
    P = transition_matrix(hmm)

    for i in 1:S
        β[T, i] = 1.0
    end
    for t in (T - 1):-1:1
        for i in 1:S
            β[t, i] = sum(P[i, j] * obs_density[t + 1, j] * β[t + 1, j] for j in 1:S) / c[t]
        end
    end
    for t in 1:T
        if all_zeros(@view β[t, :])
            println(β[t, :])
            throw(OverflowError("Log probabilities are too small in backward step."))
        end
    end
    return nothing
end

function forward_backward!(
    α::AbstractMatrix,
    β::AbstractMatrix,
    c::AbstractVector,
    γ::AbstractMatrix,
    ξ::AbstractArray{<:Real,3},
    obs_density::AbstractMatrix,
    hmm::HMM,
)
    T, S = size(obs_density, 1), nb_states(hmm)
    P = transition_matrix(hmm)

    forward!(α, c, obs_density, hmm)
    backward!(β, c, obs_density, hmm)

    for t in 1:T
        for i in 1:S
            γ[t, i] = α[t, i] * β[t, i]
        end
        sumγ = sum(@view γ[t, :])
        for i in 1:S
            γ[t, i] /= sumγ
        end
    end

    for t in 1:(T - 1)
        for i in 1:S, j in 1:S
            ξ[t, i, j] = α[t, i] * P[i, j] * obs_density[t + 1, j] * β[t + 1, j]
        end
        sumξ = sum(@view ξ[t, :, :])
        for i in 1:S, j in 1:S
            ξ[t, i, j] /= sumξ
        end
    end

    logL = sum(log.(c))  # TODO: @view

    return logL
end

## Likelihood of observations

function update_obs_density!(
    obs_density::AbstractMatrix, observations::AbstractVector, hmm::HMM
)
    T, S = length(observations), nb_states(hmm)
    for t in 1:T
        for s in 1:S
            obs_density[t, s] = densityof(emission(hmm, s), observations[t])
        end
    end
    for t in 1:T
        if all_zeros(@view obs_density[t, :])
            throw(OverflowError("Densities are too small for observations."))
        end
    end
end
