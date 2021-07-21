"""
    MultivariatePoissonProcess

A homegeneous multivariate Poisson process with integer marks.

# Fields
- `logλ::Vector{Float64}`: logarithms of the event rates.
"""
struct MultivariatePoissonProcess <: PointProcess{Int}
    logλ::Vector{Float64}
end

function default_param(::Type{MultivariatePoissonProcess}, h::History{Int})
    return ComponentVector(logλ = zeros(maximum(h.marks)))
end

function intensity(
    ::Type{MultivariatePoissonProcess},
    θ::Parameter,
    h::History{Int},
    t::Float64,
    m::Int,
)
    return exp(θ.logλ[m])
end

function mark_distribution(
    ::Type{MultivariatePoissonProcess},
    θ::Parameter,
    h::History{Int},
    t::Float64,
)
    return Categorical(exp.(θ.logλ) / sum(exp.(θ.logλ)))
end

function ground_intensity(
    ::Type{MultivariatePoissonProcess},
    θ::Parameter,
    h::History{Int},
    t::Float64,
)
    return sum(exp.(θ.logλ))
end

function ground_intensity_bound(
    ::Type{MultivariatePoissonProcess},
    θ::Parameter,
    h::History{Int},
    t::Float64,
)
    return ground_intensity(MultivariatePoissonProcess, θ, h, t), Inf
end

function integrated_ground_intensity(
    ::Type{MultivariatePoissonProcess},
    θ::Parameter,
    h::History{Int},
)
    return sum(exp.(θ.logλ)) * duration(h)
end
