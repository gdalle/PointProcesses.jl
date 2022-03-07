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
