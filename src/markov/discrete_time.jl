"""
    DiscreteMarkovChain{R<:Real}

Discrete-time Markov chain.

# Fields
- `π0::Vector{R}`: initial state distribution
- `P::Matrix{R}`: state transition matrix.
"""
@with_kw struct DiscreteMarkovChain{R<:Real} <: AbstractMarkovChain
    π0::Vector{R}
    P::Matrix{R}
end

Base.eltype(::Type{<:DiscreteMarkovChain}) = Vector{Int}

## Access

initial_distribution(dmc::DiscreteMarkovChain) = dmc.π0
transition_matrix(dmc::DiscreteMarkovChain) = dmc.P

## Simulation

function Base.rand(rng::AbstractRNG, dmc::DiscreteMarkovChain, T::Int)
    states = Vector{Int}(undef, T)
    states[1] = rand(rng, Categorical(dmc.π0))
    for t = 2:T
        states[t] = rand(rng, Categorical(dmc.P[states[t-1], :]))  # TODO: @view
    end
    return states
end

Base.rand(dmc::DiscreteMarkovChain, T::Int) = rand(Random.GLOBAL_RNG, dmc, T)

## Fitting

"""
    ContinuousMarkovChainStats

Sufficient statistics for the likelihood of a ContinuousMarkovChain.
"""
@with_kw struct DiscreteMarkovChainStats
    initialization::Vector{Float64}
    transition_count::Matrix{Float64}
end

function Distributions.suffstats(::Type{<:DiscreteMarkovChain}, states::Vector{Int})
    S = maximum(states)
    initialization = collect(1:S) .== states[1]
    transition_count = zeros(Int, S, S)
    for t = 1:length(states)-1
        transition_count[states[t], states[t+1]] += 1
    end
    return DiscreteMarkovChainStats(initialization, transition_count)
end

function Distributions.suffstats(
    ::Type{<:DiscreteMarkovChain},
    logγ::Matrix{Float64},
    logξ::Array{Float64, 3},
)
    T, S, _ = size(logξ)
    initialization = exp.(logγ[1, :])
    transition_count = zeros(Float64, S, S)
    for i = 1:S, j = 1:S
        transition_count[i, j] = sum(exp(x) for x in @view logξ[:, i, j])
    end
    return DiscreteMarkovChainStats(initialization, transition_count)
end

function Distributions.fit_mle(::Type{<:DiscreteMarkovChain}, ss::DiscreteMarkovChainStats)
    π0 = ss.initialization
    P = ss.transition_count ./ sum(ss.transition_count, dims = 2)
    return DiscreteMarkovChain(π0, P)
end

## Asymptotics

function stationary_distribution(dmc::DiscreteMarkovChain)
    π = real.(eigvecs(matrix(transition_matrix(dmc))')[:, end])
    return π / sum(π)
end
