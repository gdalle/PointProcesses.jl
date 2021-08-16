"""
    HiddenMarkovModel{Tr<:DiscreteMarkovChain,Em}

Hidden Markov Model with arbitrary transition model (must be a discrete Markov chain) and emission distributions.

# Fields
- `transitions::Tr`: state evolution process.
- `emissions::Vector{Em}`: one emission distribution per state.
"""
@with_kw struct HiddenMarkovModel{Tr<:DiscreteMarkovChain,Em}
    transitions::Tr
    emissions::Vector{Em}
end

## Access

transitions(hmm::HiddenMarkovModel) = hmm.transitions
transition_matrix(hmm::HiddenMarkovModel) = transition_matrix(transitions(hmm))
initial_distribution(hmm::HiddenMarkovModel) = initial_distribution(transitions(hmm))

emissions(hmm::HiddenMarkovModel) = hmm.emissions
emission(hmm::HiddenMarkovModel, s::Int) = hmm.emissions[s]
nstates(hmm::HiddenMarkovModel) = length(emissions(hmm))

## Simulation

function Base.rand(rng::AbstractRNG, hmm::HiddenMarkovModel, T::Int)
    states = rand(rng, transitions(hmm), T)
    observations = [rand(rng, emission(hmm, states[t])) for t = 1:T]
    return states, observations
end

Base.rand(hmm::HiddenMarkovModel, T::Int) = rand(Random.GLOBAL_RNG, hmm, T)
