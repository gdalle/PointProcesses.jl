"""
    TemporalHawkesProcess{R}

Multivariate temporal Hawkes process.

# Fields
- `λ::Vector{R}`: base event rates
- `α::Matrix{R}`: excitation amplitudes
- `β::Matrix{R}`: excitation decays
"""
@with_kw struct TemporalHawkesProcess{R} <: MultivariateTemporalPointProcess
    λ::Vector{R}
    α::Matrix{R}
    β::Matrix{R}
end

function build_transform(pp::TemporalHawkesProcess)
    M = length(pp.λ)
    return as((
        λ = as(Vector, asℝ₊, M),
        α = as(Matrix, asℝ₊, (M, M)),
        β = as(Matrix, asℝ₊, (M, M)),
    ))
end

function all_marks(pp::TemporalHawkesProcess)
    return 1:length(pp.λ)
end

function intensity(pp::TemporalHawkesProcess, h::TemporalHistory, t, m)
    error("not implemented")
end

function ground_intensity_bound(pp::TemporalHawkesProcess, h::TemporalHistory, t)
    λg = ground_intensity(pp, h, t)
    return λg, λg / 2
end
