struct ProductDistribution{D} <: AbstractMeasure
    marginals::Vector{D}
end

marginals(dist::ProductDistribution) = dist.marginals
marginal(dist::ProductDistribution, i::Integer) = dist.marginals[i]

length(dist::ProductDistribution) = length(marginals(dist))

## Simulation

function rand(rng::AbstractRNG, dist::ProductDistribution)
    return [rand(rng, marginal(dist, i)) for i = 1:length(dist)]
end

rand(dist::ProductDistribution) = rand(GLOBAL_RNG, dist)

## Likelihood

function logdensity(dist::ProductDistribution, x::AbstractVector)
    return sum(logdensity(marginal(dist, i), x[i]) for i = 1:length(dist))
end

## Prior likelihood

function logdensity(prior::ProductDistribution, dist::ProductDistribution)
    return sum(logdensity(marginal(prior, i), marginal(dist, i)) for i = 1:length(dist))
end

## Sufficient statistics

struct ProductSufficientStats{T<:Dists.SufficientStats}
    ss_marginals::Vector{T}
end

function suffstats(::Type{<:ProductDistribution{T}}, x::AbstractMatrix) where {T}
    return ProductSufficientStats([suffstats(T, x[i, :] for i in size(x, 1))])
end

function suffstats(
    ::Type{<:ProductDistribution{T}},
    x::AbstractMatrix,
    w::AbstractVector,
) where {T}
    return ProductSufficientStats([suffstats(T, x[i, :], w) for i in size(x, 1)])
end

function suffstats(
    ::Type{<:ProductDistribution{T}},
    prior::ProductDistribution{S},
    x::AbstractMatrix,
) where {T,S}
    return ProductSufficientStats([
        suffstats(T, marginal(prior, i), x[i, :]) for i in size(x, 1)
    ])
end

function suffstats(
    ::Type{<:ProductDistribution{T}},
    prior::ProductDistribution{S},
    x::AbstractMatrix,
    w::AbstractVector,
) where {T,S}
    return ProductSufficientStats([
        suffstats(T, marginal(prior, i), x[i, :], w) for i in size(x, 1)
    ])
end

## Fitting

function fit(pdtype::Type{<:ProductDistribution}, args...; kwargs...)
    return fit_mle(pdtype, args..., kwargs...)
end

function fit_mle(pdtype::Type{<:ProductDistribution}, args...; kwargs...)
    ss = suffstats(pdtype, args...; kwargs...)
    return fit_mle(pdtype, ss)
end
