"""
    PoissonProcess{R<:Real,D,M}

Homogeneous temporal Poisson process with arbitrary mark distribution.

# Fields
- `λ::R`: event rate.
- `mark_dist`: mark distribution.
"""
@with_kw struct PoissonProcess{M,R<:Real,D} <: TemporalPointProcess{M}
    λ::R
    mark_dist::D
end

intensity(pp::PoissonProcess, m=nothing) = isnothing(m) ? pp.λ : pp.λ * pdf(pp.mark_dist, m)

function Base.rand(rng::AbstractRNG, pp::PoissonProcess, tmin, tmax)
    N = rand(Poisson(pp.λ * (tmax - tmin)))
    times = rand(Uniform(tmin, tmax), N)
    marks = rand(pp.mark_dist, N)
    return TemporalHistory(times, marks, tmin, tmax)
end

function Distributions.fit(
    ::Type{PoissonProcess{M}},
    history::TemporalHistory{M},
) where {R,D,M}
    λ = nb_events(history) / duration(history)
    mark_dist = fit(D, history.marks)
    return PoissonProcess(λ, mark_dist)
end
