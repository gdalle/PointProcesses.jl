"""
    MultivariateHawkesProcess{R}

Homogeneous temporal multivariate Hawkes process.

# Fields
- `λ::Vector{R}`: base event rates
- `α::Matrix{R}`: excitation amplitudes
- `β::Matrix{R}`: excitation decays
"""
@with_kw struct MultivariateHawkesProcess{R} <: TemporalPointProcess{Int}
    λ::Vector{R}
    α::Matrix{R}
    β::Matrix{R}
end

function build_transform(pp::MultivariateHawkesProcess)
    M = length(pp.λ)
    return as((
        λ = as(Vector, asℝ₊, M),
        α = as(Matrix, asℝ₊, (M, M)),
        β = as(Matrix, asℝ₊, (M, M)),
    ))
end
