struct DiscreteProductDistribution{T<:DiscreteUnivariateDistribution} <:
        DiscreteMultivariateDistribution
    marginals::Vector{T}
end

marginals(dist::ProductDistribution) = dist.marginals
marginal(dist::ProductDistribution, i::Integer) = dist.marginals[i]

Base.length(dist::ProductDistribution) = length(marginals(dist))

## Simulation

function Distributions._rand!(
    rng::AbstractRNG,
    dist::ProductDistribution,
    x::AbstractVector,
)
    for i = 1:length(dist)
        x[i] = rand(rng, marginal(dist, i))
    end
    return x
end

function Distributions._rand!(dist::ProductDistribution, x::AbstractVector)
    return rand!(Random.GLOBAL_RNG, dist, x)
end

## Likelihood

function Distributions.logpdf(dist::ProductDistribution, x::AbstractVector)
    l = 0.0
    for i = 1:length(dist)
        l += logpdf(marginal(dist, i), x[i])
    end
    return l
end

## Prior likelihood

function Distributions.logpdf(prior::ProductDistribution, dist::ProductDistribution)
    l = 0.0
    for i in length(dist)
        l += logpdf(marginal(prior, i), marginal(dist, i))
    end
    return l
end

## Sufficient statistics

struct ProductSufficientStats{T<:SufficientStats}
    ss_marginals::Vector{T}
end

function Distributions.suffstats(
    ::Type{<:ProductDistribution{T}},
    x::AbstractMatrix,
) where {T}
    return ProductSufficientStats([suffstats(T, x[i, :] for i in size(x, 1))])
end

function Distributions.suffstats(
    ::Type{<:ProductDistribution{T}},
    x::AbstractMatrix,
    w::AbstractVector,
) where {T}
    return ProductSufficientStats([suffstats(T, x[i, :], w) for i in size(x, 1)])
end

function Distributions.suffstats(
    ::Type{<:ProductDistribution{T}},
    prior::ProductDistribution{S},
    x::AbstractMatrix,
) where {T,S}
    return ProductSufficientStats([
        suffstats(T, marginal(prior, i), x[i, :]) for i in size(x, 1)
    ])
end

function Distributions.suffstats(
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

function Distributions.fit(pdtype::Type{<:ProductDistribution}, args...; kwargs...)
    return fit_mle(pdtype, args..., kwargs...)
end

function Distributions.fit_mle(pdtype::Type{<:ProductDistribution}, args...; kwargs...)
    ss = suffstats(pdtype, args...; kwargs...)
    return fit_mle(pdtype, ss)
end
