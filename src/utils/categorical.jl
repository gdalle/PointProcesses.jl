## Warning: type piracy

## Prior likelihood

function DensityInterface.logdensityof(prior::Dirichlet, dist::Categorical)
    return logdensity(prior, probs(dist))
end

## Sufficient statistics with prior

function suffstats(
    ::Type{<:Categorical},
    prior::Dirichlet,
    x::AbstractArray{<:Integer},
)
    return CategoricalStats(add_categorical_counts!(prior.alpha .- 1, x))
end

function suffstats(
    ::Type{<:Categorical},
    prior::Dirichlet,
    x::AbstractArray{<:Integer},
    w::AbstractArray{Float64},
)
    return CategoricalStats(add_categorical_counts!(prior.alpha .- 1, x, w))
end
