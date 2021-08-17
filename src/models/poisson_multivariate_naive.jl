"""
    NaiveMultivariatePoissonProcess{R}

Homogeneous temporal multivariate Poisson process.

This naive implementation demonstrates the use of our general `TemporalPointProcess` interface.

# Fields
- `λ::Vector{R}`: event rates.
"""
struct NaiveMultivariatePoissonProcess{R} <: TemporalPointProcess{Int}
    λ::Vector{R}
end

## Conversion utilities

NaiveMultivariatePoissonProcess{R}(nt::NamedTuple) where {R} = NaiveMultivariatePoissonProcess(nt.λ)

function build_transform(pp::NaiveMultivariatePoissonProcess)
    M = length(pp.λ)
    return as((λ = as(Vector, asℝ₊, M),))
end

## Intensity functions

function intensity(pp::NaiveMultivariatePoissonProcess, h::TemporalHistory, t, m)
    return pp.λ[m]
end

function mark_distribution(pp::NaiveMultivariatePoissonProcess, h::TemporalHistory, t)
    return Categorical(pp.λ / sum(pp.λ))
end

function ground_intensity(pp::NaiveMultivariatePoissonProcess, h::TemporalHistory, t)
    return sum(pp.λ)
end

function ground_intensity_bound(pp::NaiveMultivariatePoissonProcess, h::TemporalHistory, t)
    λg = ground_intensity(pp, h, t)
    return λg, Inf
end
