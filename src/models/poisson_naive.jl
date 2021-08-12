"""
    NaivePoissonProcess{R}

Homogeneous temporal multivariate Poisson process.

This naive implementation demonstrates the use of our general `TemporalPointProcess` interface.

# Fields
- `λ::Vector{R}`: event rates.
"""
@with_kw struct NaivePoissonProcess{R} <: TemporalPointProcess{Int}
    λ::Vector{R}
end

## Conversion utilities

NaivePoissonProcess{R}(nt::NamedTuple) where {R} = NaivePoissonProcess(nt.λ)

function build_transform(pp::NaivePoissonProcess)
    M = length(pp.λ)
    return as((λ = as(Vector, asℝ₊, M),))
end

## Intensity functions

function intensity(pp::NaivePoissonProcess, h::TemporalHistory{Int}, t, m::Integer)
    return pp.λ[m]
end

function mark_distribution(pp::NaivePoissonProcess, h::TemporalHistory{Int}, t)
    return Categorical(pp.λ / sum(pp.λ))
end

function ground_intensity(pp::NaivePoissonProcess, h::TemporalHistory{Int}, t)
    return sum(pp.λ)
end

function ground_intensity_bound(pp::NaivePoissonProcess, h::TemporalHistory{Int}, t)
    λg = ground_intensity(pp, h, t)
    return λg, Inf
end
