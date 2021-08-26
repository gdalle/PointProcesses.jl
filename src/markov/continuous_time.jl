"""
    ContinuousMarkovChain{T1,T2}

Continuous-time Markov chain with finite state space.

# Fields
- `π0::T1`: initial state distribution
- `Q::T2`: rate matrix.
"""
struct ContinuousMarkovChain{T1<:AbstractVector{<:Real},T2<:AbstractMatrix{<:Real}} <:
       AbstractMarkovChain
    π0::T1
    Q::T2
end

MeasureTheory.sampletype(::ContinuousMarkovChain) = TemporalHistory{Int}

## Access

nb_states(mc::ContinuousMarkovChain) = length(mc.π0)

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
    transitions = [Dists.Categorical(P[s, :]) for s = 1:nb_states(mc)]  # TODO: @view
    waiting_times = [Dists.Exponential(1 / D[s]) for s = 1:nb_states(mc)]

    h = TemporalHistory(Float64[], Int[], tmin, tmax)
    s = rand(rng, Dists.Categorical(mc.π0))
    push!(h, tmin, s)
    t = tmin
    while t < tmax
        Δt = rand(rng, waiting_times[s])
        if t + Δt < tmax
            t += Δt
            s = rand(rng, transitions[s])
            push!(h, t, s)
        else
            break
        end
    end
    return h
end

Base.rand(mc::ContinuousMarkovChain, tmin, tmax) = rand(GLOBAL_RNG, mc, tmin, tmax)

## Prior

"""
    ContinuousMarkovChainPrior

Define a Dirichlet prior on the initial distribution and a Gamma prior on the transition rates from each state.

# Fields
- `π0α`: Dirichlet parameter
- `Pα`: Gamma rate parameters
- `Pβ`: Gamma shape parameters
"""
struct ContinuousMarkovChainPrior <: AbstractMarkovChainPrior
    π0α::Vector{Float64}
    Pα::Matrix{Float64}
    Pβ::Vector{Float64}
end

## Prior logpdf

function MeasureTheory.logdensity(
    prior::ContinuousMarkovChainPrior,
    mc::ContinuousMarkovChain,
)
    l = logdensity(Dists.Dirichlet(prior.π0α), Dists.Categorical(mc.π0))
    Q = rate_matrix(mc)
    for i = 1:nb_states(mc), j = 1:nb_states(mc)
        if i != j
            l += logdensity(Dists.Gamma(prior.Pα[i, i], 1 / prior.Pβ[i, j]), Q[i, j])
        end
    end
    return l
end

## Fitting

"""
    ContinuousMarkovChainStats

Sufficient statistics for the likelihood of a ContinuousMarkovChain.
"""
struct ContinuousMarkovChainStats <: Dists.SufficientStats
    initialization::Vector{Float64}
    transition_count::Matrix{Float64}
    duration::Vector{Float64}
end

function Dists.suffstats(::Type{ContinuousMarkovChain}, h::TemporalHistory{<:Integer})
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
    return ContinuousMarkovChainStats(initialization, transition_count, duration)
end

function Dists.suffstats(::Type{ContinuousMarkovChain}, m̂, D̂)
    error("not implemented")
end

function Dists.suffstats(
    ::Type{<:ContinuousMarkovChain},
    prior::ContinuousMarkovChainPrior,
    args...,
)
    ss = suffstats(ContinuousMarkovChain, args...)
    ss.initialization .+= (prior.π0α .- 1)
    ss.transition_count .+= (prior.Pα .- 1)
    ss.duration .+= prior.Pβ
    return ss
end

function Dists.fit_mle(::Type{<:ContinuousMarkovChain}, ss::ContinuousMarkovChainStats)
    π0 = ss.initialization
    π0 ./= sum(π0)
    Q = ss.transition_count ./ ss.duration
    Q[diagind(Q)] .= -vec(sum(Q, dims = 2)) + Q[diagind(Q)]
    return ContinuousMarkovChain(π0, Q)
end

## Asymptotics

function stationary_distribution(mc::ContinuousMarkovChain)
    π = real.(eigvecs(rate_matrix(mc)')[:, end])
    return π / sum(π)
end

## Misc

function discretize_chain(mc::ContinuousMarkovChain, Δt)
    return DiscreteMarkovChain(initial_distribution(mc), exp(rate_matrix(mc) * Δt))
end
