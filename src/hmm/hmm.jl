"""
    HMM{TransitionsType,EmissionsType}

Hidden Markov Model with arbitrary transition model (must be a discrete Markov chain) and emission distributions.

# Fields
- `transitions::TransitionsType`: state evolution process.
- `emissions::Vector{EmissionsType}`: one emission distribution per state.
"""
struct HiddenMarkovModel{TransitionsType,EmissionsType}
    transitions::TransitionsType
    emissions::Vector{EmissionsType}
end

struct HiddenMarkovModelPrior{TransitionsPriorType,EmissionsPriorType}
    transitions_prior::TransitionsPriorType
    emissions_prior::Vector{EmissionsPriorType}
end

const HMM = HiddenMarkovModel
const HMMPrior = HiddenMarkovModelPrior

## Access

transitions(hmm::HMM) = hmm.transitions
initial_distribution(hmm::HMM) = initial_distribution(transitions(hmm))
transition_matrix(hmm::HMM) = transition_matrix(transitions(hmm))

emissions(hmm::HMM) = hmm.emissions
emission(hmm::HMM, s::Integer) = hmm.emissions[s]
nb_states(hmm::HMM) = length(emissions(hmm))

## Simulation

function Base.rand(rng::AbstractRNG, hmm::HMM, T::Integer)
    states = rand(rng, transitions(hmm), T)
    observations = [rand(rng, emission(hmm, states[t])) for t = 1:T]
    return states, observations
end

Base.rand(hmm::HMM, T::Integer) = rand(GLOBAL_RNG, hmm, T)
