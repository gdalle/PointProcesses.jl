const CategoricalPrior = Dirichlet

## Prior likelihood

function Distributions.logpdf(prior::CategoricalPrior, dist::Categorical)
    return logpdf(prior, probs(dist))
end

## MAP fitting

function fit_map(
    ::Type{<:Categorical},
    prior::CategoricalPrior,
    x::AbstractVector{Integer},
)
    k = length(prior)
    counts = Distributions.add_categorical_counts!(zeros(k), x)
    corrected_counts = counts .+ prior.alpha .- 1.0
    p = Distributions.pnormalize!(corrected_counts)
    return Categorical(p)
end

function fit_map(
    ::Type{<:Categorical},
    prior::CategoricalPrior,
    x::AbstractVector{Integer},
    w::AbstractVector{Real},
)
    k = length(prior)
    counts = Distributions.add_categorical_counts!(zeros(k), x, w)
    corrected_counts = counts .+ prior.alpha .- 1.0
    p = Distributions.pnormalize!(corrected_counts)
    return Categorical(p)
end
