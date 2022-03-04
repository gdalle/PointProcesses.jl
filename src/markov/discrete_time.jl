## Structs

"""
    DiscreteMarkovChain

Discrete-time Markov chain with finite state space.

# Fields
- `π0`: initial state distribution.
- `P`: state transition matrix.
"""
Base.@kwdef struct DiscreteMarkovChain{
    R1<:Real,R2<:Real,V1<:AbstractVector{R1},M2<:AbstractMatrix{R2}
}
    π0::V1
    P::M2
end

@inline DensityInterface.DensityKind(::DiscreteMarkovChain) = HasDensity()

"""
    DiscreteMarkovChainPrior

Define a Dirichlet prior on the initial distribution and on the transitions from each state.

# Fields
- `π0α`: Dirichlet parameter for the initial distribution
- `Pα`: Dirichlet parameters for the transition matrix
"""
Base.@kwdef struct DiscreteMarkovChainPrior{
    R1<:Real,R2<:Real,V1<:AbstractVector{R1},M2<:AbstractMatrix{R2}
}
    π0α::V1
    Pα::M2
end

@inline DensityInterface.DensityKind(::DiscreteMarkovChainPrior) = HasDensity()

"""
    DiscreteMarkovChainStats

Sufficient statistics for the likelihood of a `DiscreteMarkovChain`.
"""
Base.@kwdef struct DiscreteMarkovChainStats{
    R1<:Real,R2<:Real,V1<:AbstractVector{R1},M2<:AbstractMatrix{R2}
}
    initialization::V1
    transition_count::M2
end

## Access

nb_states(mc::DiscreteMarkovChain) = length(mc.π0)
initial_distribution(mc::DiscreteMarkovChain) = mc.π0
transition_matrix(mc::DiscreteMarkovChain) = mc.P

## Simulation

function Base.rand(rng::AbstractRNG, mc::DiscreteMarkovChain, T::Integer)
    states = Vector{Int}(undef, T)
    states[1] = rand(rng, Categorical(mc.π0))
    transitions = [Categorical(mc.P[s, :]) for s in 1:nb_states(mc)]
    for t in 2:T
        states[t] = rand(rng, transitions[states[t - 1]])
    end
    return states
end

function Base.rand(rng::AbstractRNG, mc::DiscreteMarkovChainPrior)
    return error("not implemented")
end

Base.rand(mc::DiscreteMarkovChain, T::Integer) = rand(GLOBAL_RNG, mc, T)
Base.rand(prior::DiscreteMarkovChainPrior) = rand(GLOBAL_RNG, prior)

## Asymptotics

function stationary_distribution(mc::DiscreteMarkovChain)
    π = real.(eigvecs(transition_matrix(mc)')[:, end])
    return π / sum(π)
end

## Sufficient statistics

function Distributions.suffstats(
    ::Type{<:DiscreteMarkovChain}, x::AbstractVector{<:Integer}
)
    S, T = maximum(x), length(x)
    initialization = [s == x[1] for s in 1:S]
    transition_count = zeros(Int, S, S)
    for t in 1:(T - 1)
        transition_count[x[t], x[t + 1]] += 1
    end
    return DiscreteMarkovChainStats(;
        initialization=initialization, transition_count=transition_count
    )
end

# function Distributions.suffstats(
#     ::Type{<:DiscreteMarkovChain}, γ::AbstractMatrix{<:Real}, ξ::AbstractArray{<:Real,3}
# )
#     T, S, _ = size(ξ)
#     initialization = γ[1, :]
#     transition_count = zeros(Float64, S, S)
#     for i in 1:S, j in 1:S
#         transition_count[i, j] = sum(@view ξ[:, i, j])
#     end
#     return DiscreteMarkovChainStats(;
#         initialization=initialization, transition_count=transition_count
#     )
# end

## Densities

function DensityInterface.logdensityof(mc::DiscreteMarkovChain, x::AbstractVector)
    T = length(x)
    l = log(mc.π0[x[1]])
    for t in 2:T
        l += log(mc.P[x[t - 1], x[t]])
    end
    return l
end

function DensityInterface.logdensityof(
    prior::DiscreteMarkovChainPrior, mc::DiscreteMarkovChain
)
    l = logdensityof(Dirichlet(prior.π0α), Categorical(mc.π0))
    for s in 1:nb_states(mc)
        l += logdensityof(Dirichlet(@view prior.Pα[s, :]), Categorical(@view mc.P[s, :]))
    end
    return l
end

## Fit

function Distributions.fit_mle(::Type{<:DiscreteMarkovChain}, ss::DiscreteMarkovChainStats)
    π0 = ss.initialization ./ sum(ss.initialization)
    P = ss.transition_count ./ sum(ss.transition_count; dims=2)
    return DiscreteMarkovChain(; π0=π0, P=P)
end

function Distributions.fit_mle(mctype::Type{<:DiscreteMarkovChain}, args...; kwargs...)
    ss = suffstats(mctype, args...; kwargs...)
    return fit_mle(mctype, ss)
end

function fit_map(mctype::Type{<:DiscreteMarkovChain}, prior::DiscreteMarkovChainPrior, args...; kwargs...)
    ss = suffstats(mctype, args...; kwargs...)
    ss.initialization .+= (prior.π0α .- 1)
    ss.transition_count .+= (prior.Pα .- 1)
    return fit_mle(mctype, ss)
end
