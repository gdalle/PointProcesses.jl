"""
    ContinuousMarkovChain

Continuous-time Markov chain with finite state space.

# Fields
- `π0`: initial state distribution
- `Q`: rate matrix.
"""
@with_kw struct ContinuousMarkovChain{
    T1<:AbstractVector{<:Real},
    T2<:AbstractMatrix{<:Real},
} <: AbstractMarkovChain
    π0::T1
    Q::T2
end

Base.eltype(::Type{<:ContinuousMarkovChain}) = TemporalHistory{Int}

## Access

nstates(mc::ContinuousMarkovChain) = length(mc.π0)

initial_distribution(mc::ContinuousMarkovChain) = mc.π0

rate_matrix(mc::ContinuousMarkovChain) = mc.Q

rate_negdiag(mc::ContinuousMarkovChain) = -diag(rate_matrix(mc))

function transition_matrix(mc::ContinuousMarkovChain)
    P = rate_matrix(mc) ./ rate_negdiag(mc)
    P[diagind(P)] .= zero(eltype(P))
    return P
end

## Simulation

function Base.rand(rng::AbstractRNG, mc::ContinuousMarkovChain, tmin, tmax)
    D = rate_negdiag(mc)
    P = transition_matrix(mc)

    h = TemporalHistory(Float64[], Int[], tmin, tmax)
    s = rand(rng, Categorical(mc.π0))
    push!(h, tmin, s)
    t = tmin
    while t < tmax
        Δt = rand(rng, Exponential(1 / D[s]))
        if t + Δt < tmax
            t += Δt
            s = rand(rng, Categorical(P[s, :]))  # TODO: @view
            push!(h, t, s)
        else
            break
        end
    end
    return h
end

Base.rand(mc::ContinuousMarkovChain, tmin, tmax) = rand(Random.GLOBAL_RNG, mc, tmin, tmax)

## Fitting

"""
    ContinuousMarkovChainStats

Sufficient statistics for the likelihood of a ContinuousMarkovChain.
"""
@with_kw struct ContinuousMarkovChainStats{
    T1<:AbstractVector{<:Real},
    T2<:AbstractVector{<:Real},
    T3<:AbstractMatrix{<:Real},
}
    initialization::T1
    duration::T2
    transition_count::T3
end

function Distributions.suffstats(::Type{ContinuousMarkovChain}, h::TemporalHistory{<:Integer})
    states = h.marks
    S = maximum(states)
    initialization = collect(1:S) .== states[1]
    duration = zeros(Float64, S)
    transition_count = zeros(Int, S, S)
    for t = 1:length(states)-1
        duration[states[t]] += h.times[t+1] - h.times[t]
        transition_count[states[t], states[t+1]] += 1
    end
    duration[states[length(states)]] += h.tmax - h.times[length(states)]
    return ContinuousMarkovChainStats(initialization, duration, transition_count)
end

function Distributions.suffstats(::Type{ContinuousMarkovChain}; m̂, D̂)
    states = h.marks
    S = maximum(states)
    initialization = collect(1:S) .== states[1]
    duration = zeros(Float64, S)
    transition_count = zeros(Int, S, S)
    for t = 1:length(states)-1
        duration[states[t]] += h.times[t+1] - h.times[t]
        transition_count[states[t], states[t+1]] += 1
    end
    duration[states[length(states)]] += h.tmax - h.times[length(states)]
    return ContinuousMarkovChainStats(initialization, duration, transition_count)
end

function Distributions.fit_mle(
    ::Type{<:ContinuousMarkovChain},
    ss::ContinuousMarkovChainStats,
)
    π0 = ss.initialization
    Q = ss.transition_count ./ ss.duration
    Q[diagind(Q)] .= -vec(sum(Q, dims = 2)) + Q[diagind(Q)]
    return ContinuousMarkovChain(π0, Q)
end

## Asymptotics

function stationary_distribution(mc::ContinuousMarkovChain)
    π = real.(eigvecs(matrix(rate_matrix(mc))')[:, end])
    return π / sum(π)
end

## Misc

function discretize(mc::ContinuousMarkovChain, Δt)
    return DiscreteMarkovChain(initial_distribution(mc), exp(rate_matrix(mc) * Δt))
end
