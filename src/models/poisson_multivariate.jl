"""
    MultivariatePoissonProcess{R}

Homogeneous multivariate temporal Poisson process.

# Fields
- `λ::Vector{R}`: event rates.
"""
Base.@kwdef struct MultivariatePoissonProcess{R<:Real} <: TemporalPointProcess{Int,Float64}
    λ::Vector{R}
end

## Prior

Base.@kwdef struct MultivariatePoissonProcessPrior{R1<:Real,R2<:Real}
    λα::Vector{R1}
    λβ::R2
end

## Sufficient statistics

Base.@kwdef struct MultivariatePoissonProcessStats{R1<:Real,R2<:Real}
    nb_events::Vector{R1}
    duration::R2
end

## Access

Base.length(pp::MultivariatePoissonProcess) = length(pp.λ)

function mark_distribution(pp::MultivariatePoissonProcess, t=nothing, h=nothing)
    return Categorical(pp.λ / sum(pp.λ))
end

intensity(pp::MultivariatePoissonProcess, m, t=nothing, h=nothing) = pp.λ[m]
log_intensity(pp::MultivariatePoissonProcess, m, t=nothing, h=nothing) = log(pp.λ[m])

ground_intensity(pp::MultivariatePoissonProcess, t=nothing, h=nothing) = sum(pp.λ)

function ground_intensity_bound(pp::MultivariatePoissonProcess, t=nothing, h=nothing)
    return ground_intensity(pp, t, h), Inf
end

function integrated_ground_intensity(pp::MultivariatePoissonProcess, h, a, b)
    return ground_intensity(pp) * (b - a)
end

## Simulation

function Base.rand(rng::AbstractRNG, pp::MultivariatePoissonProcess, tmin::Real, tmax::Real)
    N = rand(rng, Poisson(ground_intensity(pp) * (tmax - tmin)))
    times = rand(rng, Uniform(tmin, tmax), N)
    mark_dist = mark_distribution(pp)
    marks = [rand(rng, mark_dist) for _ in 1:N]
    return History(; times=times, marks=marks, tmin=tmin, tmax=tmax)
end

function Base.rand(pp::MultivariatePoissonProcess, tmin::Real, tmax::Real)
    return rand(GLOBAL_RNG, pp, tmin, tmax)
end

## Prior likelihood

function DensityInterface.logdensityof(
    prior::MultivariatePoissonProcessPrior, pp::MultivariatePoissonProcess
)
    lλ = 0.0
    for m in 1:length(pp)
        lλ += logdensityof(Gamma(prior.λα[m], 1 / prior.λβ), pp.λ[m])
    end
    return lλ
end

## Sufficient statistics

function Distributions.suffstats(::Type{<:MultivariatePoissonProcess}, h::History)
    d = duration(h)
    ne = zeros(Int, maximum(event_marks(h)))
    for m in event_marks(h)
        ne[m] += 1
    end
    return MultivariatePoissonProcessStats(; nb_events=ne, duration=d)
end

function Distributions.suffstats(
    ::Type{<:MultivariatePoissonProcess},
    histories::AbstractVector{<:History},
    weights::AbstractVector
)
    d = 0.0
    max_mark = 0
    for h in histories
        for m in event_marks(h)
            max_mark = max(max_mark, m)
        end
    end
    ne = zeros(Float64, max_mark)
    for (h, w) in zip(histories, weights)
        d += w*duration(h)
        for m in event_marks(h)
            ne[m] += w
        end
    end
    return MultivariatePoissonProcessStats(; nb_events=ne, duration=d)
end

## Fitting

function Distributions.fit_mle(
    ::Type{<:MultivariatePoissonProcess}, ss::MultivariatePoissonProcessStats
) where {D}
    λ = ss.nb_events ./ ss.duration
    return MultivariatePoissonProcess(; λ=λ)
end

function Distributions.fit_mle(pptype::Type{<:MultivariatePoissonProcess}, args...; kwargs...)
    ss = suffstats(pptype, args...; kwargs...)
    return fit_mle(pptype, ss)
end

function fit_map(
    pptype::Type{<:MultivariatePoissonProcess},
    prior::MultivariatePoissonProcessPrior,
    args...;
    kwargs...,
)
    ss = suffstats(pptype, args...; kwargs...)
    ss.nb_events .+= prior.λα .- 1
    ss.duration += prior.λβ
    return fit_mle(pptype, ss)
end
