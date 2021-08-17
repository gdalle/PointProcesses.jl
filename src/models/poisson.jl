"""
    PoissonProcess{R<:Real,D,M}

Homogeneous temporal Poisson process with arbitrary mark distribution.

# Fields
- `λ::R`: event rate.
- `mark_dist`: mark distribution.
"""
struct PoissonProcess{M,D,R<:Real} <: TemporalPointProcess{M}
    λ::R
    mark_dist::D
end

function PoissonProcess{M}(λ::R, mark_dist::D) where {M,D,R}
    return PoissonProcess{M,D,R}(λ, mark_dist)
end

mark_distribution(pp::PoissonProcess) = pp.mark_dist

intensity(pp::PoissonProcess, m = nothing) =
    isnothing(m) ? pp.λ : pp.λ * pdf(pp.mark_dist, m)

function Base.rand(rng::AbstractRNG, pp::PoissonProcess, tmin, tmax)
    N = rand(rng, Poisson(pp.λ * (tmax - tmin)))
    times = rand(rng, Uniform(tmin, tmax), N)
    marks = rand(rng, pp.mark_dist, N)
    return TemporalHistory(times, marks, tmin, tmax)
end

function Distributions.fit(
    ::Type{<:PoissonProcess{M,D}},
    history::TemporalHistory{M},
) where {M,D}
    λ = nb_events(history) / duration(history)
    mark_dist = fit(D, history.marks)
    return PoissonProcess{M}(λ, mark_dist)
end
