abstract type AbstractMarkovChain end

"""
    DiscreteMarkovChain

Discrete-time Markov chain.

# Fields
- `π0::Vector{Float64}`: initial state distribution
- `P::Matrix{Float64}`: state transition matrix.
"""
struct DiscreteMarkovChain <: AbstractMarkovChain
    π0::Vector{Float64}
    P::Matrix{Float64}
end

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

nstates(mc::AbstractMarkovChain) = length(mc.π0)
matrix(dmc::DiscreteMarkovChain) = dmc.P
matrix(cmc::ContinuousMarkovChain) = cmc.Q

# Simulation

function Base.rand(rng::AbstractRNG, dmc::DiscreteMarkovChain, T::Int; π0)
    s = Vector{Int}(undef, T)
    s[1] = rand(Categorical(π0))
    for t = 2:T
        s[t] = rand(Categorical(@view dmc.P[s[t-1], :]))
    end
    return s
end

Base.rand(dmc::DiscreteMarkovChain, T::Int; π0 = dmc.π0) = rand(GLOBAL_RNG, dmc, T, π0 = π0)

function Base.rand(rng::AbstractRNG, cmc::ContinuousMarkovChain, tmin, tmax; π0)
    jump_rates = -diag(cmc.Q)
    jump_proba = cmc.Q ./ jump_rates
    jump_proba[diagind(jump_proba)] .= 0

    h = History{Int}(Float64[], Int[], tmin, tmax)
    s = rand(Categorical(cmc.π0))
    push!(h, t_min, s)
    t = t_min
    while t < t_max
        Δt = rand(Exponential(1 / jump_rates[s]))
        if t + Δt < t_max
            t += Δt
            s = rand(Categorical(@view jump_proba[s, :]))
            push!(h, t, s)
        else
            break
        end
    end
    return h
end

Base.rand(cmc::ContinuousMarkovChain, tmin, tmax; π0 = cmc.π0) =
    rand(GLOBAL_RNG, cmc, tmin, tmax, π0 = π0)

# Asymptotics

function stationary_distribution(mc::AbstractMarkovChain)
    π = real.(eigvecs(matrix(mc)')[:, end])
    return π / sum(π)
end

# Discretization

function discretize(cmc::ContinuousMarkovChain, Δt)
    return DiscreteMarkovChain(cmc.π0, exp(cmc.Q * Δt))
end
