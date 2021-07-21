"""
    HiddenMarkovModel{Tr, Em}

Hidden Markov Model with arbitrary transition model (of type `Tr`) and emission distributions (of type `Em`).

# Fields
- `transitions::Tr`: state evolution process.
- `emissions::Vector{Em}`: one emission distribution per state.

# Examples

```jldoctest
Random.seed!(63)
dmc = DiscreteMarkovChain([0.3, 0.7], [0.9 0.1; 0.2 0.8])
emission1 = TimedPointProcess(MultivariatePoissonProcess([0., 1., 2.]), 0., 1.)
emission2 = TimedPointProcess(MultivariatePoissonProcess([2., 1., 0.]), 0., 1.)
hmm = HiddenMarkovModel(dmc, [emission1, emission2])
states, observations = rand(hmm, 100)
sum(nb_events(observations[t]) for t = 1:100)

# output

1017
```
"""
mutable struct HiddenMarkovModel{Tr,Em}
    transitions::Tr
    emissions::Vector{Em}
end

nstates(hmm::HiddenMarkovModel) = length(hmm.emissions)

# Simulation

function Base.rand(rng::AbstractRNG, hmm::HiddenMarkovModel, T::Int)
    states = rand(rng, hmm.transitions, T)
    observations = [rand(rng, hmm.emissions[states[t]]) for t = 1:T]
    return states, observations
end

Base.rand(hmm::HiddenMarkovModel, T::Int) = rand(GLOBAL_RNG, hmm, T)

# Likelihood of observations

function update_observation_likelihood!(
    obs_logpdf::Matrix,
    hmm::HiddenMarkovModel,
    observations::Vector,
)
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
