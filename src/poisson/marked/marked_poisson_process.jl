"""
    MarkedPoissonProcess{R,M,D}

Homogeneous temporal Poisson process with arbitrary mark distribution.

# Fields

- `λ::R`: event rate.
- `mark_dist::D`: mark distribution with sample type `M`.
"""
Base.@kwdef struct MarkedPoissonProcess{R<:Real,M,D} <: AbstractPoissonProcess{M}
    λ::R
    mark_dist::D
end

function Base.show(io::IO, pp::MarkedPoissonProcess)
    return print(io, "MarkedPoissonProcess($(pp.λ), $(pp.mark_dist))")
end

## Eltype from distribution

function MarkedPoissonProcess(λ::R, mark_dist::D) where {R,D<:UnivariateDistribution}
    M = eltype(mark_dist)
    return MarkedPoissonProcess{R,M,D}(λ, mark_dist)
end

function MarkedPoissonProcess(λ::R, mark_dist::D) where {R,D<:MultivariateDistribution}
    M = Vector{eltype(mark_dist)}
    return MarkedPoissonProcess{R,M,D}(λ, mark_dist)
end

## Access

ground_intensity(pp::MarkedPoissonProcess) = pp.λ
mark_distribution(pp::MarkedPoissonProcess) = pp.mark_dist
