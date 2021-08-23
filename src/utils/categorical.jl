## Warning: type piracy

## Prior likelihood

function MeasureTheory.logdensity(prior::Dists.Dirichlet, dist::Dists.Categorical)
    return logdensity(prior, probs(dist))
end

## Sufficient statistics with prior

function Dists.suffstats(
    ::Type{<:Dists.Categorical},
    prior::Dists.Dirichlet,
    x::AbstractArray{<:Integer},
)
    return Dists.CategoricalStats(add_categorical_counts!(prior.alpha .- 1, x))
end

function Dists.suffstats(
    ::Type{<:Dists.Categorical},
    prior::Dists.Dirichlet,
    x::AbstractArray{<:Integer},
    w::AbstractArray{Float64},
)
    return Dists.CategoricalStats(add_categorical_counts!(prior.alpha .- 1, x, w))
end
