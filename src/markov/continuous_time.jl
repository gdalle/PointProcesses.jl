"""
    ContinuousMarkovChain

Continuous-time Markov chain.

# Fields
- `π0::Vector{Float64}`: initial state distribution
- `Q::Matrix{Float64}`: rate matrix.
"""
struct ContinuousMarkovChain <: AbstractMarkovChain
    π0::Vector{Float64}
    Q::Matrix{Float64}
end

Base.eltype(::Type{ContinuousMarkovChain}) = History{Int}

matrix(cmc::ContinuousMarkovChain) = cmc.Q

# Simulation

function Base.rand(
    rng::AbstractRNG,
    cmc::ContinuousMarkovChain,
    tmin::Float64,
    tmax::Float64,
)
    jump_rates = -diag(cmc.Q)
    jump_proba = cmc.Q ./ jump_rates
    jump_proba[diagind(jump_proba)] .= 0

    h = History(Float64[], Int[], tmin, tmax)
    s = rand(rng, Categorical(cmc.π0))
    push!(h, tmin, s)
    t = tmin
    while t < tmax
        Δt = rand(rng, Exponential(1 / jump_rates[s]))
        if t + Δt < tmax
            t += Δt
            s = rand(rng, Categorical(jump_proba[s, :]))  # TODO: @view
            push!(h, t, s)
        else
            break
        end
    end
    return h
end

function Base.rand(cmc::ContinuousMarkovChain, tmin::Float64, tmax::Float64)
    return rand(GLOBAL_RNG, cmc, tmin, tmax)
end

# Fitting

"""
    ContinuousMarkovChainStats

Sufficient statistics for the likelihood of a ContinuousMarkovChain.
"""
struct ContinuousMarkovChainStats
    initialization::Vector{Float64}
    duration::Vector{Float64}
    transition_count::Matrix{Float64}
end

function Distributions.suffstats(::Type{ContinuousMarkovChain}, h::History{Int})
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
    ::Type{ContinuousMarkovChain},
    ss::ContinuousMarkovChainStats,
)
    π0 = ss.initialization
    Q = ss.transition_count ./ ss.duration
    Q[diagind(Q)] .= 0.0
    Q[diagind(Q)] .= -vec(sum(Q, dims = 2))
    return ContinuousMarkovChain(π0, Q)
end

# Asymptotics

function stationary_distribution(cmc::ContinuousMarkovChain)
    π = real.(eigvecs(matrix(cmc.Q)')[:, end])
    return π / sum(π)
end

# Misc

function discretize(cmc::ContinuousMarkovChain, Δt)
    return DiscreteMarkovChain(cmc.π0, exp(cmc.Q * Δt))
end
