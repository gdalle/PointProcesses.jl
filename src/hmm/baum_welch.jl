## Likelihood of observations

function update_observation_likelihood!(
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

## Log version

function baum_welch_step!(
    logα::AbstractMatrix,
    logβ::AbstractMatrix,
    logγ::AbstractMatrix,
    logξ::AbstractArray{<:Real,3},
    obs_logpdf::AbstractMatrix,
    hmm::HiddenMarkovModel,
    observations::AbstractVector,
)
    update_observation_likelihood!(obs_logpdf, hmm, observations)
    logL = forward_backward_log!(logα, logβ, logγ, logξ, obs_logpdf, hmm)
    new_transitions = fit(typeof(transitions(hmm)), logγ, logξ)
    new_emissions =
        [fit(typeof(emission(hmm, s)), observations, exp.(logγ[:, s])) for s = 1:nstates(hmm)]  # TODO: @view
    new_hmm = HiddenMarkovModel(new_transitions, new_emissions)
    return new_hmm, logL
end

function baum_welch!(
    logα::AbstractMatrix,
    logβ::AbstractMatrix,
    logγ::AbstractMatrix,
    logξ::AbstractArray{<:Real,3},
    obs_logpdf::AbstractMatrix,
    hmm::HiddenMarkovModel,
    observations::AbstractVector;
    iter,
)
    logL_evolution = Float64[]
    for _ = 1:iter
        hmm, logL = baum_welch_step!(logα, logβ, logγ, logξ, obs_logpdf, hmm, observations)
        push!(logL_evolution, logL)
    end
    return hmm, logL_evolution
end

## Non mutating version

function baum_welch(hmm::HiddenMarkovModel, observations::AbstractVector; iter, log = true)
    T, S = length(observations), nstates(hmm)
    if log
        logα = Matrix{Float64}(undef, T, S)
        logβ = Matrix{Float64}(undef, T, S)
        logγ = Matrix{Float64}(undef, T, S)
        logξ = Array{Float64,3}(undef, T-1, S, S)
        obs_logpdf = Matrix{Float64}(undef, T, S)
        hmm_est, logL_evolution =
            baum_welch!(logα, logβ, logγ, logξ, obs_logpdf, hmm, observations; iter = iter)
    else
        error("Not implemented yet")
    end
    return hmm_est, logL_evolution
end
