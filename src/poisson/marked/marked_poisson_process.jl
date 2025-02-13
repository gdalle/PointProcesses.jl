"""
    MarkedPoissonProcess{R,D}

Homogeneous temporal Poisson process with arbitrary mark distribution.

# Fields

- `λ::R`: event rate.
- `mark_dist::D`: mark distribution.

# Constructor

    MarkedPoissonProcess(λ, mark_dist)
"""
struct MarkedPoissonProcess{R<:Real,D} <: AbstractPoissonProcess
    λ::R
    mark_dist::D
end

function Base.show(io::IO, pp::MarkedPoissonProcess)
    return print(io, "MarkedPoissonProcess($(pp.λ), $(pp.mark_dist))")
end

## Access

ground_intensity(pp::MarkedPoissonProcess) = pp.λ
mark_distribution(pp::MarkedPoissonProcess) = pp.mark_dist
