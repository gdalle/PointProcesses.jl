"""
    DiscreteMarkovChain

Discrete-time Markov chain..

# Fields
- `π0::Vector{Float64}`: initial state distribution
- `P::Matrix{Float64}`: state transition matrix.

# Examples

```jldoctest
using Random; Random.seed!(63)
dmc = DiscreteMarkovChain([0.3, 0.7], [0.9 0.1; 0.2 0.8])
states = rand(dmc, 100)
fit(DiscreteMarkovChain, states)

# output

DiscreteMarkovChain([1.0, 0.0], [0.9425287356321839 0.05747126436781609; 0.4166666666666667 0.5833333333333334])
```
"""
struct DiscreteMarkovChain <: AbstractMarkovChain
    π0::Vector{Float64}
    P::Matrix{Float64}
end

Base.eltype(::Type{DiscreteMarkovChain}) = Vector{Int}

# Simulation

function Base.rand(rng::AbstractRNG, dmc::DiscreteMarkovChain, T::Int)
    states = Vector{Int}(undef, T)
    states[1] = rand(rng, Categorical(dmc.π0))
    for t = 2:T
        states[t] = rand(rng, Categorical(dmc.P[states[t-1], :]))  # TODO: @view
    end
    return states
end

function Base.rand(dmc::DiscreteMarkovChain, T::Int)
    return rand(GLOBAL_RNG, dmc, T)
end

# Fitting

"""
    ContinuousMarkovChainStats

Sufficient statistics for the likelihood of a ContinuousMarkovChain.
"""
struct DiscreteMarkovChainStats
    initialization::Vector{Float64}
    transition_count::Matrix{Float64}
end

function Distributions.suffstats(::Type{DiscreteMarkovChain}, states::Vector{Int})
    S = maximum(states)
    initialization = collect(1:S) .== states[1]
    transition_count = zeros(Int, S, S)
    for t = 1:length(states)-1
        transition_count[states[t], states[t+1]] += 1
    end
    return DiscreteMarkovChainStats(initialization, transition_count)
end

function Distributions.fit_mle(::Type{DiscreteMarkovChain}, ss::DiscreteMarkovChainStats)
    π0 = ss.initialization
    P = ss.transition_count ./ sum(ss.transition_count, dims = 2)
    return DiscreteMarkovChain(π0, P)
end

# Asymptotics

function stationary_distribution(dmc::DiscreteMarkovChain)
    π = real.(eigvecs(matrix(dmc.P)')[:, end])
    return π / sum(π)
end
