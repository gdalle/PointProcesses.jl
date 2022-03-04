"""
    PoissonProcess{D,M,R}

Homogeneous temporal Poisson process with arbitrary mark distribution.

# Fields
- `λ::R`: event rate.
- `mark_dist::D`: mark distribution with sample type `M`.
"""
Base.@kwdef struct PoissonProcess{D,M,R<:Real} <: TemporalPointProcess{M,Float64}
    λ::R
    mark_dist::D
end

function PoissonProcess(λ::R, mark_dist::D) where {R,D<:UnivariateDistribution}
    M = eltype(mark_dist)
    return PoissonProcess{D,M,R}(λ, mark_dist)
end

function PoissonProcess(λ::R, mark_dist::D) where {R,D<:MultivariateDistribution}
    M = Vector{eltype(mark_dist)}
    return PoissonProcess{D,M,R}(λ, mark_dist)
end

## Prior

Base.@kwdef struct PoissonProcessPrior{D,R1<:Real,R2<:Real}
    λα::R1
    λβ::R2
    mark_dist_prior::D
end

## Sufficient statistics

Base.@kwdef struct PoissonProcessStats{R1<:Real,R2<:Integer,MDSS}
    duration::R1
    nb_events::R2
    mark_dist_suffstats::MDSS
end

## Access

mark_distribution(pp::PoissonProcess) = pp.mark_dist
ground_intensity(pp::PoissonProcess) = pp.λ

function intensity(pp::PoissonProcess, m)
    return ground_intensity(pp) * densityof(mark_distribution(pp), m)
end

function log_intensity(pp::PoissonProcess, m)
    return log(ground_intensity(pp)) + logdensityof(mark_distribution(pp), m)
end

## Simulation

function Base.rand(rng::AbstractRNG, pp::PoissonProcess, tmin::Real, tmax::Real)
    N = rand(rng, Poisson(ground_intensity(pp) * (tmax - tmin)))
    times = rand(rng, Uniform(tmin, tmax), N)
    marks = [rand(rng, mark_distribution(pp)) for _ in 1:N]
    return History(; times=times, marks=marks, tmin=tmin, tmax=tmax)
end

Base.rand(pp::PoissonProcess, tmin::Real, tmax::Real) = rand(GLOBAL_RNG, pp, tmin, tmax)

## Likelihood

function DensityInterface.logdensityof(pp::PoissonProcess, h::History)
    l = -ground_intensity(pp) * duration(h)
    mark_dist = mark_distribution(pp)
    for m in event_marks(h)
        l += logdensityof(mark_dist, m)
    end
    return l
end

## Prior likelihood

function DensityInterface.logdensityof(prior::PoissonProcessPrior, pp::PoissonProcess)
    lλ = logdensityof(Gamma(prior.λα, 1 / prior.λβ), pp.λ)
    ld = logdensityof(prior.mark_dist_prior, pp.mark_dist)
    return lλ + ld
end

## Sufficient statistics

function Distributions.suffstats(::Type{<:PoissonProcess{D}}, h::History) where {D}
    return PoissonProcessStats(;
        duration=duration(h),
        nb_events=nb_events(h),
        mark_dist_suffstats=suffstats(D, event_marks(h)),
    )
end

## Fitting

function Distributions.fit_mle(
    ::Type{<:PoissonProcess{D}}, ss::PoissonProcessStats
) where {D}
    λ = ss.nb_events / ss.duration
    mark_dist = fit_mle(D, ss.mark_dist_suffstats)
    return PoissonProcess(; λ=λ, mark_dist=mark_dist)
end

function Distributions.fit_mle(pptype::Type{<:PoissonProcess}, args...; kwargs...)
    ss = suffstats(pptype, args...; kwargs...)
    return fit_mle(pptype, ss)
end

function fit_map(
    ::Type{<:PoissonProcess{D}}, prior::PoissonProcessPrior, h::History
) where {D}
    λ = (nb_events(h) + prior.λα - 1) / (duration(h) + prior.λβ)
    mark_dist = fit_map(D, prior.mark_dist_prior, history.marks)
    return PoissonProcess(λ, mark_dist)
end
