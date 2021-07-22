abstract type MultivariatePointProcess <: PointProcess{Int} end

function all_marks(pptype::Type{<:MultivariatePointProcess}, θ::Parameter) end

all_marks(pp::MultivariatePointProcess) = all_marks(typeof(pp), params(pp))

function ground_intensity(
    pptype::Type{<:MultivariatePointProcess},
    θ::Parameter,
    h::History{Int},
    t::Float64,
)
    return sum(intensity(pptype, θ, h, t, m) for m in all_marks(pptype, θ))
end

function all_mark_probabilities(
    pptype::Type{<:MultivariatePointProcess},
    θ::Parameter,
    h::History{Int},
    t::Float64,
)
    λg = ground_intensity(pptype, θ, h, t)
    return [intensity(pptype, θ, h, t, m) / λg for m in all_marks(pptype, θ)]
end

function mark_distribution(
    pptype::Type{<:MultivariatePointProcess},
    θ::Parameter,
    h::History{Int},
    t::Float64,
)
    return DiscreteNonParametric(
        all_marks(pptype, θ),
        all_mark_probabilities(pptype, θ, h, t),
    )
end
