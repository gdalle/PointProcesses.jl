"""
    MarkedPoissonProcess{M,R,D}

Homogeneous temporal Poisson process with arbitrary mark distribution.

# Fields

- `λ::R`: event rate.
- `mark_dist::D`: mark distribution with sample type `M`.

# Constructor

    MarkedPoissonProcess{M}(λ, mark_dist)
"""
Base.@kwdef struct MarkedPoissonProcess{M,R<:Real,D} <: AbstractPoissonProcess{M}
    λ::R
    mark_dist::D

    function MarkedPoissonProcess{M}(λ::R, mark_dist::D) where {M,R,D}
        return new{M,R,D}(λ, mark_dist)
    end
end

function Base.show(io::IO, pp::MarkedPoissonProcess)
    return print(io, "MarkedPoissonProcess($(pp.λ), $(pp.mark_dist))")
end

## Eltype from distribution

"""
    MarkedPoissonProcess(λ, mark_dist::UnivariateDistribution)

Construct a `MarkedPoissonProcess` with scalar marks.
"""
function MarkedPoissonProcess(λ, mark_dist::UnivariateDistribution)
    M = eltype(mark_dist)
    return MarkedPoissonProcess{M}(λ, mark_dist)
end

"""
    MarkedPoissonProcess(λ, mark_dist::MultivariateDistribution)

Construct a `MarkedPoissonProcess` with vector marks.
"""
function MarkedPoissonProcess(λ, mark_dist::MultivariateDistribution)
    M = Vector{eltype(mark_dist)}
    return MarkedPoissonProcess{M}(λ, mark_dist)
end

## Access

ground_intensity(pp::MarkedPoissonProcess) = pp.λ
mark_distribution(pp::MarkedPoissonProcess) = pp.mark_dist
