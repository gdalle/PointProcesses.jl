"""
    MvCategorical

Product of independent `Categorical` distributions.
"""
struct MvCategorical{R<:Real,T<:AbstractVector{<:AbstractVector{R}}} <:
       DiscreteMultivariateDistribution
    p::T
end

## Attributes

Base.length(dist::MvCategorical) = length(dist.p)

## Simulation

function Distributions._rand!(
    rng::AbstractRNG,
    dist::MvCategorical,
    x::AbstractVector{<:Real},
)
    for i = 1:length(dist)
        x[i] = rand(rng, Categorical(dist.p[i]))
    end
    return x
end

## Likelihood

function Distributions._logpdf(dist::MvCategorical, x::AbstractVector{Integer})
    l = 0.0
    for i = 1:length(dist)
        l += log(dist.p[i][x[i]])
    end
    return l
end

## Prior

struct MvCategoricalPrior{R<:Real,T<:AbstractVector{<:AbstractVector{R}}}
    α::T
end

function Distributions.logpdf(prior::MvCategoricalPrior, dist::MvCategorical)
    l = 0.0
    for i in length(dist)
        l += logpdf(CategoricalProir(prior.α[i]), Categorical(dist.p[i]))
    end
    return l
end

## MLE

function Distributions.fit_mle(
    ::Type{<:MvCategorical},
    ks::AbstractVector{Integer},
    x::AbstractMatrix{Integer},
)
    marginals = [fit_mle(Categorical, ks[i], x[i, :]) for i = 1:length(ks)]
    p = [probs(dist) for dist in marginals]
    return MvCategorical(p)
end

function Distributions.fit_mle(
    ::Type{<:MvCategorical},
    ks::AbstractVector{Integer},
    x::AbstractMatrix{Integer},
    w::AbstractVector{Real},
)
    marginals = [fit_mle(Categorical, ks[i], x[i, :], w) for i = 1:length(ks)]
    p = [probs(dist) for dist in marginals]
    return MvCategorical(p)
end


function Distributions.fit_mle(::Type{<:MvCategorical}, x::AbstractMatrix{Integer})
    ks = vec(maximum(x, dims = 2))
    return fit_mle(MvCategorical, ks, x)
end

function Distributions.fit_mle(
    ::Type{<:MvCategorical},
    x::AbstractMatrix{Integer},
    w::AbstractVector{Real},
)
    ks = vec(maximum(x, dims = 2))
    return fit_mle(MvCategorical, ks, x, w)
end

## MAP estimation

function fit_map(
    ::Type{<:MvCategorical},
    prior::MvCategoricalPrior,
    x::AbstractMatrix{Integer},
)
    marginals = [
        fit_map(Categorical, CategoricalPrior(prior.α[i]), x[i, :]) for
        i = 1:length(prior.α)
    ]
    p = [probs(dist) for dist in marginals]
    return MvCategorical(p)
end

function fit_map(
    ::Type{<:MvCategorical},
    prior::MvCategoricalPrior,
    x::AbstractMatrix{Integer},
    w::AbstractVector{Real},
)
    marginals = [
        fit_map(Categorical, CategoricalPrior(prior.α[i]), x[i, :], w) for
        i = 1:length(prior.α)
    ]
    p = [probs(dist) for dist in marginals]
    return MvCategorical(p)
end
