"""
    HMM{E}

Hidden Markov Model with arbitrary emission distributions of type `E`.

# Fields
- `transitions::DiscreteMarkovChain`: state evolution process.
- `emissions::Vector{E}`: one emission distribution per state.
"""
struct HMM{E}
    transitions::DiscreteMarkovChain
    emissions::Vector{E}
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

function observations_logpdf(hmm::HMM, observations)
    T, S = length(observations), nstates(hmm)
    obs_logpdf = Matrix{Float64}(undef, T, S)
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
    return obs_logpdf
end
