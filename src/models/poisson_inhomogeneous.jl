"""
    InhomogeneousPoissonProcess{D,M}

Inhomogeneous temporal Poisson process with arbitrary mark distribution.

# Fields
- `λ::Function`: intensity function.
- `mark_dist::D`: mark distribution.
"""
@with_kw struct InhomogeneousPoissonProcess{D,M} <: TemporalPointProcess{M}
    λ::Function
    mark_dist::D
end
