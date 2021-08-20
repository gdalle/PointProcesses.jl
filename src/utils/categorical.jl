## Prior likelihood

function Distributions.logpdf(prior::Dirichlet, dist::Categorical)
    return logpdf(prior, probs(dist))
end

## Sufficient statistics with prior

function Distributions.suffstats(::Type{<:Categorical}, prior::Dirichlet, x::AbstractArray{<:Integer})
    return CategoricalStats(add_categorical_counts!(prior.alpha .- 1, x))
end

function Distributions.suffstats(::Type{<:Categorical}, prior::Dirichlet, x::AbstractArray{<:Integer}, w::AbstractArray{Float64})
    return CategoricalStats(add_categorical_counts!(prior.alpha .- 1, x, w))
end
