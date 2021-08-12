"""
    PoissonProcess{R<:Real,D,M}

Homogeneous temporal Poisson process with arbitrary mark distribution.

# Fields
- `λ::R`: event rate.
- `mark_dist`: mark distribution.
"""
@with_kw struct PoissonProcess{R<:Real,D,M} <: TemporalPointProcess{M}
    λ::R
    mark_dist::D
end

function Base.rand(rng::AbstractRNG, pp::PoissonProcess, tmin, tmax)
    N = rand(Poisson(pp.λ * (tmax - tmin)))
    times = rand(Uniform(tmin, tmax), N)
    marks = rand(pp.mark_dist, N)
    return TemporalHistory(times, marks, tmin, tmax)
end

function Distributions.fit(
    ::Type{PoissonProcess{R,D,M}},
    history::TemporalHistory,
) where {R,D,M}
    λ = nb_events(history) / duration(history)
    mark_dist = fit(D, history.marks)
    return PoissonProcess(λ, mark_dist)
end
