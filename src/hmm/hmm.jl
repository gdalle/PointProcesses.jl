"""
    HMM{Tr<:DiscreteMarkovChain,Em}

Hidden Markov Model with arbitrary transition model (must be a discrete Markov chain) and emission distributions.

# Fields
- `transitions::Tr`: state evolution process.
- `emissions::Vector{Em}`: one emission distribution per state.
"""
struct HiddenMarkovModel{Tr<:DiscreteMarkovChain,Em} <: AbstractMeasure
    transitions::Tr
    emissions::Vector{Em}
end

const HMM = HiddenMarkovModel

## Access

transitions(hmm::HMM) = hmm.transitions
initial_distribution(hmm::HMM) = initial_distribution(transitions(hmm))
transition_matrix(hmm::HMM) = transition_matrix(transitions(hmm))

emissions(hmm::HMM) = hmm.emissions
emission(hmm::HMM, s::Int) = hmm.emissions[s]
nb_states(hmm::HMM) = length(emissions(hmm))

## Simulation

function Base.rand(rng::AbstractRNG, hmm::HMM, T::Integer)
    states = rand(rng, transitions(hmm), T)
    observations = [rand(rng, emission(hmm, states[t])) for t = 1:T]
    return states, observations
end

Base.rand(hmm::HMM, T::Integer) = rand(GLOBAL_RNG, hmm, T)
