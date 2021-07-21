"""
    HMM{Tr, Em, Ob}

Hidden Markov Model with arbitrary transition model (of type `Tr`), emission distributions (of type `Em`) and observations (of type `Ob`).

# Fields
- `transitions::Tr`: state evolution process.
- `emissions::Vector{Em}`: one emission distribution per state with element type `Ob`.
"""
mutable struct HMM{Tr, Em, Ob}
    transitions::Tr
    emissions::Vector{Em}
end

nstates(hmm::HMM) = length(hmm.emissions)

# Simulation

function Base.rand(rng::AbstractRNG, hmm::HMM, T)
    states = rand(hmm.transitions, T)
    observations = [rand(hmm.emissions[states[t]]) for t = 1:T]
    return states, observations
end

Base.rand(hmm::HMM, T) = rand(GLOBAL_RNG, hmm, T)

# Likelihood of observations

function update_observation_likelihood!(obs_logpdf::Matrix, hmm::HMM, observations::Vector)
    T, S = length(observations), nstates(hmm)
    for t = 1:T
        for s = 1:S
            obs_logpdf[t, s] = logpdf(hmm.emissions[s], observations[t])
        end
    end
    for t = 1:T
        if all_minus_inf(@view obs_logpdf[t, :])
            throw(OverflowError("Log-probabilities are too small for observations."))
        end
    end
end
