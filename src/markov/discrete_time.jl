"""
    DiscreteMarkovChain{T1,T2}

Discrete-time Markov chain with finite state space.

# Fields
- `π0::T1`: initial state distribution.
- `P::T2`: state transition matrix.
"""
struct DiscreteMarkovChain{T1<:AbstractVector{<:Real},T2<:AbstractMatrix{<:Real}} <:
       AbstractMarkovChain
    π0::T1
    P::T2
end

MeasureTheory.sampletype(::DiscreteMarkovChain) = Vector{Int}

## Access

nb_states(mc::DiscreteMarkovChain) = length(mc.π0)

initial_distribution(mc::DiscreteMarkovChain) = mc.π0

transition_matrix(mc::DiscreteMarkovChain) = mc.P

## Simulation

function Base.rand(rng::AbstractRNG, mc::DiscreteMarkovChain, T::Integer)
    states = Vector{Int}(undef, T)
    states[1] = rand(rng, Dists.Categorical(mc.π0))
    transitions = [Dists.Categorical(mc.P[s, :]) for s = 1:nb_states(mc)]  # TODO: @view
    for t = 2:T
        states[t] = rand(rng, transitions[states[t-1]])
    end
    return states
end

Base.rand(mc::DiscreteMarkovChain, T::Integer) = rand(GLOBAL_RNG, mc, T)

## Logpdf

function MeasureTheory.logdensity(mc::DiscreteMarkovChain, states::Vector{Integer})
    T = length(states)
    l = log(mc.π0[states[1]])
    for t = 2:T
        l += log(mc.P[states[t-1], states[t]])
    end
    return l
end

## Prior

"""
    DiscreteMarkovChainPrior

Define a Dirichlet prior on the initial distribution and on the transitions from each state.

# Fields
- `π0α`: Dirichlet parameter for the initial distribution
- `Pα`: Dirichlet parameters for the transition matrix
"""
struct DiscreteMarkovChainPrior <: AbstractMarkovChainPrior
    π0α::Vector{Float64}
    Pα::Matrix{Float64}
end

## Prior logpdf

function MeasureTheory.logdensity(prior::DiscreteMarkovChainPrior, mc::DiscreteMarkovChain)
    l = logdensity(Dists.Dirichlet(prior.π0α), Dists.Categorical(mc.π0))
    for s = 1:nb_states(mc)
        l += logdensity(Dists.Dirichlet(prior.Pα[s, :]), Dists.Categorical(mc.P[s, :]))
    end
    return l
end

## Fitting

"""
    DiscreteMarkovChainStats

Sufficient statistics for the likelihood of a `DiscreteMarkovChain`.
"""
struct DiscreteMarkovChainStats <: Dists.SufficientStats
    initialization::Vector{Float64}
    transition_count::Matrix{Float64}
end

function Dists.suffstats(::Type{<:DiscreteMarkovChain}, states::Vector{<:Integer})
    S = maximum(states)
    initialization = collect(1:S) .== states[1]
    transition_count = zeros(Int, S, S)
    for t = 1:length(states)-1
        transition_count[states[t], states[t+1]] += 1
    end
    return DiscreteMarkovChainStats(initialization, transition_count)
end

function Dists.suffstats(
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

function Dists.suffstats(
    ::Type{<:DiscreteMarkovChain},
    prior::DiscreteMarkovChainPrior,
    args...,
)
    ss = suffstats(DiscreteMarkovChain, args...)
    ss.initialization .+= (prior.π0α .- 1)
    ss.transition_count .+= (prior.Pα .- 1)
    return ss
end

function Dists.fit_mle(::Type{<:DiscreteMarkovChain}, ss::DiscreteMarkovChainStats)
    π0 = ss.initialization
    π0 ./= sum(π0)
    P = ss.transition_count
    P ./= sum(ss.transition_count, dims = 2)
    return DiscreteMarkovChain(π0, P)
end

## Asymptotics

function stationary_distribution(mc::DiscreteMarkovChain)
    π = real.(eigvecs(transition_matrix(mc)')[:, end])
    return π / sum(π)
end
