function forward_nolog!(
    α::AbstractMatrix,
    c::AbstractVector,
    obs_pdf::AbstractMatrix,
    hmm::HiddenMarkovModel,
)
    T, S = size(obs_logpdf, 1), nstates(hmm)
    π0, P = initial_distribution(hmm), transition_matrix(hmm)
    for i = 1:S
        α[1, i] = π0[i] * obs_pdf[1, i]
    end
    c[1] = 1 / sum(@view α[1, :])  # scaling
    for i = 1:S
        α[1, i] *= c[1]
    end
    for t = 1:T-1
        for j = 1:S
            α[t+1, j] = sum(α[t, i] * P[i, j]) * obs_pdf[t+1, j]
        end
        c[t+1] = 1 / sum(@view α[t+1, :])  # scaling
        for j = 1:S
            α[t+1, j] *= c[t+1]
        end
    end
    for t = 1:T
        if all_zeros(@view α[t, :])
            throw(OverflowError("Probabilities are too small in forward step."))
        end
    end
    return logα
end

function forward_log!(
    logα::AbstractMatrix,
    obs_logpdf::AbstractMatrix,
    hmm::HiddenMarkovModel,
)
    T, S = size(obs_logpdf, 1), nstates(hmm)
    logπ0, logP = log.(initial_distribution(hmm)),
    log.(transition_matrix(hmm))
    for i = 1:S
        logα[1, i] = logπ0[i] + obs_logpdf[1, i]
    end
    for t = 1:T-1
        for j = 1:S
            logα[t+1, j] = logsumexp(logα[t, :] + logP[:, j]) + obs_logpdf[t+1, j]  # TODO
        end
    end
    for t = 1:T
        if all_nan(@view logα[t, :])
            throw(OverflowError("Log probabilities are too small in forward step."))
        end
    end
    return logα
end

function backward_nolog!(
    β::AbstractMatrix,
    c::AbstractVector,
    obs_pdf::AbstractMatrix,
    hmm::HiddenMarkovModel,
)
    T, S = size(obs_logpdf, 1), nstates(hmm)
    P = transition_matrix(hmm)

    for i = 1:S
        β[T, i] = 0.0
    end
    for t = T-1:-1:1
        for i = 1:S
            β[t, i] = sum(P[i, j] * obs_pdf[t+1, j] * β[t+1, j] for j = 1:S) * c[t]
        end
    end
    for t = 1:T
        if all_zeros(@view β[t, :])
            throw(OverflowError("Log probabilities are too small in backward step."))
        end
    end
end


function backward_log!(
    logβ::AbstractMatrix,
    obs_logpdf::AbstractMatrix,
    hmm::HiddenMarkovModel,
)
    T, S = size(obs_logpdf, 1), nstates(hmm)
    logP = log.(transition_matrix(hmm))

    for i = 1:S
        logβ[T, i] = 0.0
    end
    for t = T-1:-1:1
        for i = 1:S
            logβ[t, i] = logsumexp(logP[i, :] .+ obs_logpdf[t+1, :] .+ logβ[t+1, :])  # TODO
        end
    end
    for t = 1:T
        if all_nan(@view logβ[t, :])
            throw(OverflowError("Log probabilities are too small in backward step."))
        end
    end
end

function forward_backward_nolog!(
    α::AbstractMatrix,
    β::AbstractMatrix,
    c::AbstractVector,
    γ::AbstractMatrix,
    ξ::AbstractArray{<:Real, 3},
    obs_pdf::AbstractMatrix,
    hmm::HiddenMarkovModel,
)
    T, S = size(obs_logpdf, 1), nstates(hmm)

    forward_nolog!(α, c, hmm, obs_pdf)
    backward_nolog!(β, c, hmm, obs_pdf)

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
            ξ[t, i, j] = α[t, i] * P[i, j] * obs_pdf[t+1, j] * β[t+1, j]
        end
        sumξ = sum(@view ξ[t, :, :])
        for i = 1:S, j = 1:S
            ξ[t, i, j] /= sumξ
        end
    end

    L = sum(@view α[T, :])

    return log(L)
end


function forward_backward_log!(
    logα::AbstractMatrix,
    logβ::AbstractMatrix,
    logγ::AbstractMatrix,
    logξ::AbstractArray{<:Real, 3},
    obs_logpdf::AbstractMatrix,
    hmm::HiddenMarkovModel,
)
    T, S = size(obs_logpdf, 1), nstates(hmm.transitions)
    logP = log.(transition_matrix(hmm))

    forward_log!(logα, obs_logpdf, hmm)
    backward_log!(logβ, obs_logpdf, hmm)

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
            logξ[t, i, j] = logα[t, i] + logP[i, j] + obs_logpdf[t+1, j] + logβ[t+1, j]
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

function update_obs_pdf!(
    obs_pdf::AbstractMatrix,
    hmm::HiddenMarkovModel,
    observations::AbstractVector,
)
    T, S = length(observations), nstates(hmm)
    for t = 1:T
        for s = 1:S
            obs_pdf[t, s] = pdf(emission(hmm, s), observations[t])
        end
    end
    for t = 1:T
        if all_zeros(@view obs_logpdf[t, :])
            throw(OverflowError("Probabilities are too small for observations."))
        end
    end
end

function update_obs_logpdf!(
    obs_logpdf::AbstractMatrix,
    hmm::HiddenMarkovModel,
    observations::AbstractVector,
)
    T, S = length(observations), nstates(hmm)
    for t = 1:T
        for s = 1:S
            obs_logpdf[t, s] = logpdf(emission(hmm, s), observations[t])
        end
    end
    for t = 1:T
        if all_minus_inf(@view obs_logpdf[t, :])
            throw(OverflowError("Log-probabilities are too small for observations."))
        end
    end
end

function baum_welch_nolog!(
    α::AbstractMatrix,
    β::AbstractMatrix,
    c::AbstractVector,
    γ::AbstractMatrix,
    ξ::AbstractArray{<:Real,3},
    obs_pdf::AbstractMatrix,
    hmm::HiddenMarkovModel,
    observations::AbstractVector;
    iterations,
)
    logL_evolution = Float64[]
    for _ = 1:iterations
        update_obs_pdf!(obs_pdf, hmm, observations)
        logL = forward_backward_nolog!(α, β, c, γ, ξ, obs_pdf, hmm)
        push!(logL_evolution, logL)
        new_transitions = fit(typeof(transitions(hmm)), γ = γ, ξ = ξ)
        new_emissions =
            [fit(typeof(emission(hmm, s)), observations, γ[:, s]) for s = 1:nstates(hmm)]  # TODO: @view
        hmm = HiddenMarkovModel(new_transitions, new_emissions)
    end
    return hmm, logL_evolution
end

function baum_welch_log!(
    logα::AbstractMatrix,
    logβ::AbstractMatrix,
    logγ::AbstractMatrix,
    logξ::AbstractArray{<:Real,3},
    obs_logpdf::AbstractMatrix,
    hmm::HiddenMarkovModel,
    observations::AbstractVector;
    iterations,
)
    logL_evolution = Float64[]
    for _ = 1:iterations
        update_obs_logpdf!(obs_logpdf, hmm, observations)
        logL = forward_backward_log!(logα, logβ, logγ, logξ, obs_logpdf, hmm)
        push!(logL_evolution, logL)
        new_transitions = fit(typeof(transitions(hmm)), logγ = logγ, logξ = logξ)
        new_emissions = [
            fit(typeof(emission(hmm, s)), observations, exp.(logγ[:, s])) for
            s = 1:nstates(hmm)
        ]  # TODO: @view
        hmm = HiddenMarkovModel(new_transitions, new_emissions)
    end
    return hmm, logL_evolution
end

## Non mutating version

function baum_welch(
    hmm::HiddenMarkovModel,
    observations::AbstractVector;
    iterations,
    log = true,
)
    T, S = length(observations), nstates(hmm)
    if log
        logα = Matrix{Float64}(undef, T, S)
        logβ = Matrix{Float64}(undef, T, S)
        logγ = Matrix{Float64}(undef, T, S)
        logξ = Array{Float64,3}(undef, T - 1, S, S)
        obs_logpdf = Matrix{Float64}(undef, T, S)
        hmm_est, logL_evolution = baum_welch!(
            logα,
            logβ,
            logγ,
            logξ,
            obs_logpdf,
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
        obs_pdf = Matrix{Float64}(undef, T, S)
        hmm_est, logL_evolution =
            baum_welch!(α, β, c, γ, ξ, obs_pdf, hmm, observations; iterations = iterations)
    end
    return hmm_est, logL_evolution
end
