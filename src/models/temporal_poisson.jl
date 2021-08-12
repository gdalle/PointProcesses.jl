"""
    TemporalPoissonProcess{R}

Homogeneous temporal multivariate Poisson process.

# Fields
- `λ::Vector{R}`: event rates.
"""
@with_kw struct TemporalPoissonProcess{R} <: MultivariateTemporalPointProcess
    λ::Vector{R}
end

TemporalPoissonProcess{R}(nt::NamedTuple) where {R} = TemporalPoissonProcess(nt.λ)

function build_transform(pp::TemporalPoissonProcess)
    M = length(pp.λ)
    return as((λ = as(Vector, asℝ₊, M),))
end

function all_marks(pp::TemporalPoissonProcess)
    return 1:length(pp.λ)
end

function intensity(pp::TemporalPoissonProcess, h::TemporalHistory{Int}, t, m::Integer)
    return pp.λ[m]
end

function ground_intensity_bound(pp::TemporalPoissonProcess, h::TemporalHistory{Int}, t)
    λg = ground_intensity(pp, h, t)
    return λg, Inf
end
