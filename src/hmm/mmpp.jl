"""
    MarkovModulatedPoissonProcess{M,Tr<:ContinuousMarkovChain,Em<:TemporalPointProcess{M}}

Markov-Modulated Poisson Process with mark type `M`.

# Fields
- `transitions::Tr`: state evolution process.
- `emissions::Vector{Em}`: one emission distribution per state.
"""
struct MarkovModulatedPoissonProcess{
    M,
    Tr<:ContinuousMarkovChain,
    Em<:TemporalPointProcess{M},
} <: AbstractMeasure
    transitions::Tr
    emissions::Vector{Em}
end

const MMPP = MarkovModulatedPoissonProcess

## Access

transitions(mmpp::MMPP) = mmpp.transitions
initial_distribution(mmpp::MMPP) = initial_distribution(transitions(mmpp))
transition_matrix(mmpp::MMPP) = transition_matrix(transitions(mmpp))
rate_matrix(mmpp::MMPP) = rate_matrix(transitions(mmpp))

emissions(mmpp::MMPP) = mmpp.emissions
emission(mmpp::MMPP, s::Int) = mmpp.emissions[s]
nb_states(mmpp::MMPP) = length(emissions(mmpp))

## Simulation

function Base.rand(rng::AbstractRNG, mmpp::MMPP{M}, tmin, tmax) where {M}
    state_history = rand(rng, transitions(mmpp), tmin, tmax)
    transition_times, states = event_times(state_history), event_marks(state_history)
    observations = TemporalHistory(Float64[], M[], tmin, tmax)
    for k = 1:length(transition_times)
        local_s = states[k]
        local_tmin = transition_times[k]
        local_tmax = k < length(transition_times) ? transition_times[k+1] : tmax
        local_observations = rand(rng, emission(mmpp, local_s), local_tmin, local_tmax)
        append!(observations, local_observations)
    end
    return state_history, observations
end

Base.rand(mmpp::MMPP, args...) = rand(GLOBAL_RNG, mmpp, args...)
