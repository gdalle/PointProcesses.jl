abstract type AbstractMarkovChain end

struct DiscreteMarkovChain <: AbstractMarkovChain
    π0::Vector{Float64}
    P::Matrix{Float64}
end

struct ContinuousMarkovChain <: AbstractMarkovChain
    π0::Vector{Float64}
    Q::Matrix{Float64}
end

nstates(mc::AbstractMarkovChain) = length(mc.π0)
matrix(dmc::DiscreteMarkovChain) = dmc.P
matrix(cmc::ContinuousMarkovChain) = cmc.Q

# Simulation

function Distributions.rand(dmc::DiscreteMarkovChain, T::Int; π0 = dmc.π0)
    s = Vector{Int}(undef, T)
    s[1] = rand(Categorical(π0))
    for t = 2:T
        s[t] = rand(Categorical(@view dmc.P[s[t-1], :]))
    end
    return s
end

function Distributions.rand(cmc::ContinuousMarkovChain, tmin, tmax; π0 = cmc.π0)
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

# Asymptotics

function stationary_distribution(mc::AbstractMarkovChain)
    π = real.(eigvecs(matrix(mc)')[:, end])
    return π / sum(π)
end

# Discretization

function discretize(cmc::ContinuousMarkovChain, Δt)
    return DiscreteMarkovChain(cmc.π0, exp(cmc.Q * Δt))
end
