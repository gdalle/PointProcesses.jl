"""
    TemporalPointProcess{L<:Real, M}

The common supertype of all temporal point processes (i.e. point processes on the real line) with mark type `M`.
"""
abstract type TemporalPointProcess{M} <: AbstractPointProcess{Float64,M} end

"""
    BoundedTemporalPointProcess{L<:Real, M}

Store a temporal point process `P` with pre-defined start and end times.
"""
@with_kw struct BoundedTemporalPointProcess{M,P<:TemporalPointProcess{M}} <: TemporalPointProcess{M}
    pp::P
    tmin::Float64
    tmax::Float64
end

#=
"""
    MultivariateTemporalPointProcess

The common supertype of all temporal point processes with integer marks.
"""
abstract type MultivariateTemporalPointProcess <: TemporalPointProcess{Int} end

"""
    all_marks(pp)

Enumerate possible marks for a multivariate temporal point process.
"""
function all_marks(pp::MultivariateTemporalPointProcess)
    error("not implemented")
end

function ground_intensity(pp::MultivariateTemporalPointProcess, h::TemporalHistory{Int}, t)
    return sum(intensity(pp, h, t, m) for m in all_marks(pp))
end

function all_mark_probabilities(
    pp::MultivariateTemporalPointProcess,
    h::TemporalHistory{Int},
    t,
)
    λg = ground_intensity(pp, h, t)
    return [intensity(pp, h, t, m) / λg for m in all_marks(pp)]
end

function mark_distribution(pp::MultivariateTemporalPointProcess, h::TemporalHistory{Int}, t)
    return DiscreteNonParametric(all_marks(pp), all_mark_probabilities(pp, h, t))
end
=#
