"""
    PoissonProcess{D,M}

Homogeneous temporal Poisson process with arbitrary mark distribution.

# Fields
- `λ::Float64`: event rate.
- `mark_dist::D`: mark distribution.
"""
struct PoissonProcess{D,M} <: TemporalPointProcess{M}
    λ::Float64
    mark_dist::D
end

function PoissonProcess(λ::Float64, mark_dist::D) where {D}
    return PoissonProcess{D,sampletype(mark_dist)}(λ, mark_dist)
end

## Access

function mark_distribution(pp::PoissonProcess)
    return pp.mark_dist
end

function intensity(pp::PoissonProcess, m = nothing)
    if isnothing(m)
        return pp.λ
    else
        return pp.λ * pdf(mark_distribution(pp), m)
    end
end

## Simulation

function Base.rand(rng::AbstractRNG, pp::PoissonProcess, tmin::Real, tmax::Real)
    N = rand(rng, Dists.Poisson(intensity(pp) * (tmax - tmin)))
    times = rand(rng, Dists.Uniform(tmin, tmax), N)
    marks = [rand(rng, mark_distribution(pp)) for _ = 1:N]
    return TemporalHistory(times, marks, tmin, tmax)
end

Base.rand(pp::PoissonProcess, tmin::Real, tmax::Real) = rand(GLOBAL_RNG, pp, tmin, tmax)

## Likelihood

function MeasureTheory.logdensity(pp::PoissonProcess, h::TemporalHistory)
    l = -intensity(pp) * duration(h)
    mark_dist = mark_distribution(pp)
    for m in event_marks(h)
        l += logdensity(mark_dist, m)
    end
    return m
end

## Prior

struct PoissonProcessPrior{D}
    λα::Float64
    λβ::Float64
    mark_dist_prior::D
end

## Prior likelihood

function MeasureTheory.logdensity(prior::PoissonProcessPrior, pp::PoissonProcess)
    lλ = logdensity(Dists.Gamma(prior.λα, 1 / prior.λβ), pp.λ)
    ld = logdensity(prior.mark_dist_prior, pp.mark_dist)
    return lλ + ld
end

## Fitting

function Dists.fit(::Type{<:PoissonProcess{D}}, history::TemporalHistory) where {D}
    λ = nb_events(history) / duration(history)
    mark_dist = fit(D, history.marks)
    return PoissonProcess(λ, mark_dist)
end

function Dists.fit(
    ::Type{<:PoissonProcess{D}},
    prior::PoissonProcessPrior,
    history::TemporalHistory,
) where {D}
    λ = (nb_events(history) + prior.α - 1) / (duration(history) + prior.β)
    mark_dist = fit(D, prior.mark_dist_prior, history.marks)
    return PoissonProcess(λ, mark_dist)
end
