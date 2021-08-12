"""
    InhomogeneousPoissonProcess{R}

Inhomogeneous temporal multivariate Poisson process.

# Fields
- `λ::Vector{R}`: event rates.
"""
@with_kw struct InhomogeneousPoissonProcess{R} <: TemporalPointProcess{Int}
    λ::Vector{R}
end
