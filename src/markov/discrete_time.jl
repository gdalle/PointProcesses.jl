"""
    DiscreteMarkovChain

Discrete-time Markov chain with finite state space.

# Fields
- `π0::AbstractVector{<:Real}`: initial state distribution
- `P::AbstractMatrix{<:Real}`: state transition matrix.
"""
@with_kw struct DiscreteMarkovChain{
    R1<:Real,
    T1<:AbstractVector{R1},
    R2<:Real,
    T2<:AbstractVector{R2},
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

function Base.rand(rng::AbstractRNG, mc::DiscreteMarkovChain, T::Int)
    states = Vector{Int}(undef, T)
    states[1] = rand(rng, Categorical(mc.π0))
    for t = 2:T
        states[t] = rand(rng, Categorical(mc.P[states[t-1], :]))  # TODO: @view
    end
    return states
end

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
    T1<:AbstractVector{Float64},
    T2<:AbstractMatrix{Float64},
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
    T2<:AbstractVector{<:Real},
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
    ::Type{<:DiscreteMarkovChain};
    γ = nothing,
    ξ = nothing,
    logγ = nothing,
    logξ = nothing,
)
    if isnothing(γ) || isnothing(ξ)
        γ, ξ = exp.(logγ), exp.(logξ)
    end
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
