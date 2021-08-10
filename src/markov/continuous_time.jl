"""
    ContinuousMarkovChain{R<:Real}

Continuous-time Markov chain.

# Fields
- `π0::Vector{R}`: initial state distribution
- `Q::Matrix{R}`: rate matrix.
"""
@with_kw struct ContinuousMarkovChain{R<:Real} <: AbstractMarkovChain
    π0::Vector{R}
    Q::Matrix{R}
end

Base.eltype(::Type{<:ContinuousMarkovChain}) = TemporalHistory{Int}

## Access

initial_distribution(cmc::ContinuousMarkovChain) = cmc.π0

rate_matrix(cmc::ContinuousMarkovChain) = cmc.Q

rate_diag(cmc::ContinuousMarkovChain) = -diag(rate_matrix(cmc))

function transition_matrix(cmc::ContinuousMarkovChain)
    P = rate_matrix(cmc) ./ rate_diag(cmc)
    P[diagind(P)] .= zero(eltype(P))
    return P
end

## Simulation

function Base.rand(rng::AbstractRNG, cmc::ContinuousMarkovChain, tmin, tmax)
    D = rate_diag(cmc)
    P = transition_matrix(cmc)

    h = TemporalHistory(Float64[], Int[], tmin, tmax)
    s = rand(rng, Categorical(cmc.π0))
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

Base.rand(cmc::ContinuousMarkovChain, tmin, tmax) = rand(Random.GLOBAL_RNG, cmc, tmin, tmax)

## Fitting

"""
    ContinuousMarkovChainStats

Sufficient statistics for the likelihood of a ContinuousMarkovChain.
"""
@with_kw struct ContinuousMarkovChainStats
    initialization::Vector{Float64}
    duration::Vector{Float64}
    transition_count::Matrix{Float64}
end

function Distributions.suffstats(::Type{ContinuousMarkovChain}, h::TemporalHistory{Int})
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
    Q[diagind(Q)] .= 0.0
    Q[diagind(Q)] .= -vec(sum(Q, dims = 2))
    return ContinuousMarkovChain(π0, Q)
end

## Asymptotics

function stationary_distribution(cmc::ContinuousMarkovChain)
    π = real.(eigvecs(matrix(rate_matrix(cmc))')[:, end])
    return π / sum(π)
end

## Misc

function discretize(cmc::ContinuousMarkovChain, Δt)
    return DiscreteMarkovChain(initial_distribution(cmc), exp(rate_matrix(cmc) * Δt))
end
