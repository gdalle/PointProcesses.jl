## Structs

"""
    ContinuousMarkovChain

Continuous-time Markov chain with finite state space.

# Fields
- `π0`: initial state distribution
- `Q`: rate matrix.
"""
Base.@kwdef struct ContinuousMarkovChain{
    R1<:Real,R2<:Real,V1<:AbstractVector{R1},M2<:AbstractMatrix{R2}
}
    π0::V1
    Q::M2
end

@inline DensityInterface.DensityKind(::ContinuousMarkovChain) = HasDensity()

"""
    ContinuousMarkovChainPrior

Define a Dirichlet prior on the initial distribution and a Gamma prior on the transition rates from each state.

# Fields
- `π0α`: Dirichlet parameter
- `Pα`: Gamma rate parameters
- `Pβ`: Gamma shape parameters
"""
Base.@kwdef struct ContinuousMarkovChainPrior{
    R1<:Real,
    R2<:Real,
    R3<:Real,
    V1<:AbstractVector{R1},
    M2<:AbstractMatrix{R2},
    M3<:AbstractMatrix{R3},
}
    π0α::V1
    Pα::M2
    Pβ::M3
end

@inline DensityInterface.DensityKind(::ContinuousMarkovChainPrior) = HasDensity()

"""
    ContinuousMarkovChainStats

Sufficient statistics for the likelihood of a ContinuousMarkovChain.
"""
Base.@kwdef struct ContinuousMarkovChainStats{
    R1<:Real,
    R2<:Real,
    R3<:Real,
    V1<:AbstractVector{R1},
    V2<:AbstractVector{R2},
    M3<:AbstractMatrix{R3},
}
    initialization::V1
    duration::V2
    transition_count::M3
end

## Access

nb_states(mc::ContinuousMarkovChain) = length(mc.π0)
initial_distribution(mc::ContinuousMarkovChain) = mc.π0
rate_matrix(mc::ContinuousMarkovChain) = mc.Q
rate_negdiag(mc::ContinuousMarkovChain) = -diag(rate_matrix(mc))

function embedded_transition_matrix(mc::ContinuousMarkovChain)
    P = rate_matrix(mc) ./ rate_negdiag(mc)
    P[diagind(P)] .= zero(eltype(P))
    return P
end

## Simulation

function Base.rand(rng::AbstractRNG, mc::ContinuousMarkovChain, tmin, tmax)
    D = rate_negdiag(mc)
    P = embedded_transition_matrix(mc)
    transitions = [Categorical(P[s, :]) for s in 1:nb_states(mc)]
    waiting_times = [Exponential(1 / D[s]) for s in 1:nb_states(mc)]

    h = History(; times=Float64[], marks=Int[], tmin=tmin, tmax=tmax)
    s = rand(rng, Categorical(mc.π0))
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

function Base.rand(rng::AbstractRNG, mc::ContinuousMarkovChainPrior)
    return error("not implemented")
end

Base.rand(mc::ContinuousMarkovChain, tmin, tmax) = rand(GLOBAL_RNG, mc, tmin, tmax)
Base.rand(prior::ContinuousMarkovChainPrior) = rand(GLOBAL_RNG, prior)

## Asymptotics

function stationary_distribution(mc::ContinuousMarkovChain)
    π = real.(eigvecs(rate_matrix(mc)')[:, end])
    return π / sum(π)
end

## Sufficient statistics

function Distributions.suffstats(::Type{ContinuousMarkovChain}, h::History{<:Integer})
    states = h.marks
    S = maximum(states)
    initialization = collect(1:S) .== states[1]
    duration = zeros(Float64, S)
    transition_count = zeros(Int, S, S)
    for t in 1:(length(states) - 1)
        duration[states[t]] += h.times[t + 1] - h.times[t]
        transition_count[states[t], states[t + 1]] += 1
    end
    duration[states[length(states)]] += h.tmax - h.times[length(states)]
    return ContinuousMarkovChainStats(;
        initialization=initialization, duration=duration, transition_count=transition_count
    )
end

## Densities

function DensityInterface.logdensityof(mc::ContinuousMarkovChain, h::History{<:Integer})
    return error("not implemented")
end

function DensityInterface.logdensityof(
    prior::ContinuousMarkovChainPrior, mc::ContinuousMarkovChain
)
    l = logdensity(Dirichlet(prior.π0α), Categorical(mc.π0))
    Q = rate_matrix(mc)
    for i in 1:nb_states(mc), j in 1:nb_states(mc)
        if i != j
            l += logdensity(Gamma(prior.Pα[i, i], 1 / prior.Pβ[i, j]), Q[i, j])
        end
    end
    return l
end

## Fit

function Distributions.fit_mle(
    ::Type{<:ContinuousMarkovChain}, ss::ContinuousMarkovChainStats
)
    π0 = ss.initialization ./ sum(ss.initialization)
    Q = ss.transition_count ./ ss.duration
    Q[diagind(Q)] .= -vec(sum(Q; dims=2)) + Q[diagind(Q)]
    return ContinuousMarkovChain(; π0=π0, Q=Q)
end

function Distributions.fit_mle(mctype::Type{<:ContinuousMarkovChain}, args...; kwargs...)
    ss = suffstats(mctype, args...; kwargs...)
    return fit_mle(mctype, ss)
end

function fit_map(mctype::Type{<:ContinuousMarkovChain}, prior::ContinuousMarkovChainPrior, args...; kwargs...)
    ss = suffstats(mctype, args...; kwargs...)
    ss.initialization .+= (prior.π0α .- 1)
    ss.transition_count .+= (prior.Pα .- 1)
    ss.duration .+= prior.Pβ
    return fit_mle(mctype, ss)
end

## Misc

function discretize_chain(mc::ContinuousMarkovChain, Δt)
    return DiscreteMarkovChain(; π0=initial_distribution(mc), P=exp(rate_matrix(mc) * Δt))
end
