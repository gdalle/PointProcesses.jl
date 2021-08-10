"""
    TemporalHawkesProcess{R}

Multivariate temporal Hawkes process.

# Fields
- `λ::Vector{R}`: base event rates
- `α::Matrix{R}`: excitation amplitudes
- `β::Matrix{R}`: excitation decays
"""
struct HawkesProcess{R} <: MultivariateTemporalPointProcess
    λ::Vector{R}
    α::Matrix{R}
    β::Matrix{R}
end

PoissonProcess{R}(nt::NamedTuple) where {R} = PoissonProcess(nt.λ, nt.α, nt.β)

function build_transform(pp::HawkesProcess)
    M = length(pp.λ)
    return as((
        λ = as(Vector, asℝ₊, M),
        α = as(Matrix, asℝ₊, (M, M)),
        β = as(Matrix, asℝ₊, (M, M)),
    ))
end

function all_marks(pp::HawkesProcess)
    return 1:length(pp.λ)
end

function intensity(pp::HawkesProcess, h::History{Int}, t::Float64, m::Int)
    error("not implemented")
end

function ground_intensity_bound(pp::HawkesProcess, h::History{Int}, t::Float64)
    λg = ground_intensity(pp, h, t)
    return λg, λg/2
end
