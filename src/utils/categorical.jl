struct CategoricalPrior
    α::Vector{Float64}
end

## Prior likelihood

function Distributions.logpdf(prior::CategoricalPrior, dist::Categorical)
    α, p = prior.α, probs(dist)
    return logpdf(Dirichlet(α), p) # TODO
    # return logpdf(Dirichlet(α[α.>1.0]), p[α.>1.0])
end

## MAP fitting

function fit_map(
    ::Type{<:Categorical},
    prior::CategoricalPrior,
    k::Integer,
    x::AbstractVector{Integer},
)
    counts = Distributions.add_categorical_counts!(zeros(k), x)
    corrected_counts = counts .+ prior.α .- 1.0
    p = Distributions.pnormalize!(corrected_counts)
    return Categorical(p)
end

function fit_map(
    ::Type{<:Categorical},
    prior::CategoricalPrior,
    k::Integer,
    x::AbstractVector{Integer},
    w::AbstractVector{Real},
)
    counts = Distributions.add_categorical_counts!(zeros(k), x, w)
    corrected_counts = counts .+ prior.α .- 1.0
    p = Distributions.pnormalize!(corrected_counts)
    return Categorical(p)
end
