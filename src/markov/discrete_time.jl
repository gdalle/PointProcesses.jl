"""
    DiscreteMarkovChain

Discrete-time Markov chain with finite state space.

# Fields
- `π0`: initial state distribution
- `P`: state transition matrix.
"""
@with_kw struct DiscreteMarkovChain{
    T1<:AbstractVector{<:Real},
    T2<:AbstractMatrix{<:Real},
} <: AbstractMarkovChain
    π0::T1
    P::T2
end

Base.eltype(::Type{<:DiscreteMarkovChain}) = Vector{Int}

## Access

nstates(mc::DiscreteMarkovChain) = length(mc.π0)

initial_distribution(mc::DiscreteMarkovChain) = mc.π0

transition_matrix(mc::DiscreteMarkovChain) = mc.P

## Simulation

function Base.rand(rng::AbstractRNG, mc::DiscreteMarkovChain, T::Integer)
    states = Vector{Int}(undef, T)
    states[1] = rand(rng, Categorical(mc.π0))
    for t = 2:T
        states[t] = rand(rng, Categorical(mc.P[states[t-1], :]))  # TODO: @view
    end
    return states
end

Base.rand(mc::DiscreteMarkovChain, T::Integer) = rand(Random.GLOBAL_RNG, mc, T)

## Logpdf

function Distributions.logpdf(mc::DiscreteMarkovChain, states::Vector{Integer})
    T = length(states)
    l = log(mc.π0[states[1]])
    for t = 2:T
        l += log(mc.P[states[t-1], states[t]])
    end
    return l
end

## Prior

@with_kw struct DiscreteMarkovChainPrior{
    T1<:AbstractVector{<:Real},
    T2<:AbstractMatrix{<:Real},
} <: AbstractMarkovChainPrior
    π0α::T1
    Pα::T2
end

## Prior logpdf

function Distributions.logpdf(prior::DiscreteMarkovChainPrior, mc::DiscreteMarkovChain)
    l = logpdf(CategoricalPrior(prior.π0α), Categorical(mc.π0))
    for s = 1:nstates(mc)
        l += logpdf(CategoricalPrior(prior.Pα[s, :]), Categorical(mc.P[s, :]))
    end
    return l
end

## Fitting

"""
    DiscreteMarkovChainStats

Sufficient statistics for the likelihood of a DiscreteMarkovChain.
"""
@with_kw struct DiscreteMarkovChainStats{
    T1<:AbstractVector{<:Real},
    T2<:AbstractMatrix{<:Real},
}
    initialization::T1
    transition_count::T2
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
    γ::AbstractMatrix,
    ξ::AbstractArray{<:Real,3},
)
    T, S, _ = size(ξ)
    initialization = γ[1, :]
    transition_count = zeros(Float64, S, S)
    for i = 1:S, j = 1:S
        transition_count[i, j] = sum(@view ξ[:, i, j])
    end
    return DiscreteMarkovChainStats(initialization, transition_count)
end

function Distributions.fit_mle(::Type{<:DiscreteMarkovChain}, ss::DiscreteMarkovChainStats)
    π0 = ss.initialization
    P = ss.transition_count ./ sum(ss.transition_count, dims = 2)
    return DiscreteMarkovChain(π0, P)
end

## Asymptotics

function stationary_distribution(mc::DiscreteMarkovChain)
    π = real.(eigvecs(matrix(transition_matrix(mc))')[:, end])
    return π / sum(π)
end
