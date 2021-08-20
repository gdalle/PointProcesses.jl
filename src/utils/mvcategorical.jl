"""
    MvCategorical

Product of independent `Categorical` distributions.
"""
const MvCategorical = ProductDistribution{Categorical}

"""
    MvCategoricalPrior

Conjugate multivariate Dirichlet prior over the parameters of a `MvCategorical` distribution.
"""
const MvCategoricalPrior = ProductDistribution{Dirichlet}

# ## MLE

# function Distributions.fit_mle(
#     ::Type{<:MvCategorical},
#     ks::AbstractVector{Integer},
#     x::AbstractMatrix{Integer},
# )
#     marginals = [fit_mle(Categorical, ks[i], x[i, :]) for i = 1:length(ks)]
#     p = [probs(dist) for dist in marginals]
#     return MvCategorical(p)
# end

# function Distributions.fit_mle(
#     ::Type{<:MvCategorical},
#     ks::AbstractVector{Integer},
#     x::AbstractMatrix{Integer},
#     w::AbstractVector{Real},
# )
#     marginals = [fit_mle(Categorical, ks[i], x[i, :], w) for i = 1:length(ks)]
#     p = [probs(dist) for dist in marginals]
#     return MvCategorical(p)
# end


# function Distributions.fit_mle(::Type{<:MvCategorical}, x::AbstractMatrix{Integer})
#     ks = vec(maximum(x, dims = 2))
#     return fit_mle(MvCategorical, ks, x)
# end

# function Distributions.fit_mle(
#     ::Type{<:MvCategorical},
#     x::AbstractMatrix{Integer},
#     w::AbstractVector{Real},
# )
#     ks = vec(maximum(x, dims = 2))
#     return fit_mle(MvCategorical, ks, x, w)
# end

# ## MAP estimation

# function fit_map(
#     ::Type{<:MvCategorical},
#     prior::MvCategoricalPrior,
#     x::AbstractMatrix{Integer},
# )
#     marginals = [
#         fit_map(Categorical, CategoricalPrior(prior.α[i]), x[i, :]) for
#         i = 1:length(prior.α)
#     ]
#     p = [probs(dist) for dist in marginals]
#     return MvCategorical(p)
# end

# function fit_map(
#     ::Type{<:MvCategorical},
#     prior::MvCategoricalPrior,
#     x::AbstractMatrix{Integer},
#     w::AbstractVector{Real},
# )
#     marginals = [
#         fit_map(Categorical, CategoricalPrior(prior.α[i]), x[i, :], w) for
#         i = 1:length(prior.α)
#     ]
#     p = [probs(dist) for dist in marginals]
#     return MvCategorical(p)
# end
