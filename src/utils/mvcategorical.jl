function uniformprobvec(n::Int)::Vector{Float64}
    return ones(n) ./ n
end

function randprobvec(n::Int)::Vector{Float64}
    p = rand(n)
    return p ./ sum(p)
end

struct MvCategorical <: DiscreteMultivariateDistribution
    marginals::Vector{Categorical{Float64, Vector{Float64}}}
end

## Attributes

function Distributions.ndims(dist::MvCategorical)::Int
    return length(dist.marginals)
end

Base.length(dist::MvCategorical) = ndims(dist)
Base.eltype(::MvCategorical) = Int

function Distributions.ncategories(dist::MvCategorical, dim::Int)::Int
    return ncategories(dist.marginals[dim])
end

function Distributions.ncategories(dist::MvCategorical)::Vector{Int}
    return [ncategories(dist, dim) for dim = 1:length(dist)]
end

function Distributions.probs(dist::MvCategorical, dim::Int)::Vector{Float64}
    return probs(dist.marginals[dim])
end

function Distributions.probs(dist::MvCategorical)::Vector{Vector{Float64}}
    return [probs(dist, dim) for dim = 1:length(dist)]
end

function prob(dist::MvCategorical, dim::Int, val::Int)::Float64
    return probs(dist, dim)[val]
end

## Check support

function Distributions.insupport(dist::MvCategorical, x::Vector{Int})::Bool
    length(x) == length(dist) || return false
    for dim = 1:length(dist)
        if !insupport(dist.marginals[dim], x[dim])
            return false
        end
    end
    return true
end

## Likelihood

function Distributions._logpdf(dist::MvCategorical, x::Vector{Int})::Float64
    insupport(dist, x) || return -Inf
    s = 0.0
    for dim = 1:length(dist)
        s += log(prob(dist, dim, x[dim]))
    end
    return s
end

## Simulation

function Distributions._rand!(
    rng::AbstractRNG,
    dist::MvCategorical,
    x::AbstractVector{T},
) where {T<:Real}
    length(x) == length(dist) || throw(DimensionMismatch("Invalid argument dimension."))
    for dim = 1:length(dist)
        x[dim] = rand(rng, Categorical(probs(dist, dim)))
    end
    return x
end

## MLE estimation

function Distributions.fit_mle(
    ::Type{<:MvCategorical},
    ranges::Vector{Int},
    x::Matrix{Int},
)::MvCategorical
    marginals = [fit_mle(Categorical, ranges[dim], x[dim, :]) for dim = 1:length(ranges)]
    return MvCategorical(marginals)
end

function Distributions.fit_mle(
    ::Type{<:MvCategorical},
    ranges::Vector{Int},
    x::Matrix{Int},
    w::Vector{Float64},
)::MvCategorical
    marginals = [fit_mle(Categorical, ranges[dim], x[dim, :], w) for dim = 1:length(ranges)]
    return MvCategorical(marginals)
end

function Distributions.fit_mle(::Type{<:MvCategorical}, x::Matrix{Int})
    @assert size(x, 2) > 0 "Zero samples provided"
    return fit_mle(MvCategorical, maximum(x, dims = 2)[:, 1], x)
end

function Distributions.fit_mle(::Type{<:MvCategorical}, x::Matrix{Int}, w::Vector{Float64})
    @assert size(x, 2) > 0 "Zero samples provided"
    return fit_mle(MvCategorical, maximum(x, dims = 2)[:, 1], x, w)
end

## Default estimation

Distributions.fit(::Type{<:MvCategorical}, x) = fit_mle(MvCategorical, x)
Distributions.fit(::Type{<:MvCategorical}, x, w) = fit_mle(MvCategorical, x, w)

## Priors

struct CategoricalPrior
    α::Vector{Float64}
end

Base.length(prior::CategoricalPrior)::Int = length(prior.α)

struct MvCategoricalPrior
    marginals::Vector{CategoricalPrior}
end

## MAP estimation

function fit_map(::Type{<:Categorical}, k::Int, x::Vector{Int}, prior::CategoricalPrior)
    Categorical(
        Distributions.pnormalize!(
            Distributions.add_categorical_counts!(zeros(k), x) .+ prior.α .- 1.0,
        ),
        check_args = false,
    )
end

function fit_map(
    ::Type{<:Categorical},
    k::Int,
    x::Vector{Int},
    w::Vector{Float64},
    prior::CategoricalPrior,
)
    Categorical(
        Distributions.pnormalize!(
            Distributions.add_categorical_counts!(zeros(k), x, w) .+ prior.α .- 1.0,
        ),
        check_args = false,
    )
end

function fit_map(
    ::Type{<:MvCategorical},
    ranges::Vector{Int},
    x::Matrix{Int},
    prior::MvCategoricalPrior,
)::MvCategorical
    marginals = [fit_map(Categorical, ranges[dim], x[dim, :], prior.marginals[dim]) for dim = 1:length(ranges)]
    return MvCategorical(marginals)
end

function fit_map(
    ::Type{<:MvCategorical},
    ranges::Vector{Int},
    x::Matrix{Int},
    w::Vector{Float64},
    prior::MvCategoricalPrior,
)::MvCategorical
    marginals = [fit_map(Categorical, ranges[dim], x[dim, :], w, prior.marginals[dim]) for dim = 1:length(ranges)]
    return MvCategorical(marginals)
end

function fit_map(::Type{<:MvCategorical}, x::Matrix{Int}, prior::MvCategoricalPrior)::MvCategorical
    @assert size(x, 2) > 0 "Zero samples provided"
    return fit_map(MvCategorical, length.(prior.marginals), x, prior)
end

function fit_map(
    ::Type{<:MvCategorical},
    x::Matrix{Int},
    w::Vector{Float64},
    α::Vector{Vector{Float64}},
)::MvCategorical
    @assert size(x, 2) > 0 "Zero samples provided"
    return fit_map(MvCategorical, length.(α), x, w, α)
end

## Prior likelihood

function Distributions.logpdf(prior::CategoricalPrior, dist::Categorical)::Float64
    p, α = probs(dist), prior.α
    # return logpdf(Dirichlet(α), p) # TODO
    return logpdf(Dirichlet(α[α.>1.0]), p[α.>1.0])
end

function Distributions.logpdf(prior::MvCategoricalPrior, dist::MvCategorical)::Float64
    LL = 0.0
    for dim = 1:length(dist)
        LL += logpdf(prior.marginals[dim], dist.marginals[dim])
    end
    return LL
end

## Plotting

# function plot_histogram(dist::MvCategorical; logscale = false)
#     fig, ax = subplots(1, ndims(d), figsize = (ndims(d), 2), sharey = true)
#     for dim = 1:ndims(d)
#         a = ndims(d) > 1 ? ax[dim] : ax
#         if logscale
#             a.set_yscale("log")
#         end
#         a.bar(1:ncategories(d, dim), probs(d, dim))
#         a.set_xlabel("Dim $dim")
#         a.set_ylim(0, 1)
#         a.set_xticks(1:ncategories(d, dim))
#         dim == 1 ? a.set_ylabel("Probability") : nothing
#     end
#     fig.suptitle("Multivariate categorical distribution")
#     fig.tight_layout()
#     savefig("p.svg")
#     show()
# end
