"""
    MultivariateTemporalPointProcess

The common supertype of all point processes on the real line with integer marks.
"""
abstract type MultivariateTemporalPointProcess <: TemporalPointProcess{Int} end

"""
    all_marks(pp)

Enumerate possible marks for a multivariate temporal point process.
"""
function all_marks(pp::MultivariateTemporalPointProcess)
    error("not implemented")
end

function ground_intensity(pp::MultivariateTemporalPointProcess, h::History{Int}, t::Float64)
    return sum(intensity(pp, h, t, m) for m in all_marks(pp))
end

function all_mark_probabilities(
    pp::MultivariateTemporalPointProcess,
    h::History{Int},
    t::Float64,
)
    λg = ground_intensity(pp, h, t)
    return [intensity(pp, h, t, m) / λg for m in all_marks(pp)]
end

function mark_distribution(
    pp::MultivariateTemporalPointProcess,
    h::History{Int},
    t::Float64,
)
    return DiscreteNonParametric(all_marks(pp), all_mark_probabilities(pp, h, t))
end
