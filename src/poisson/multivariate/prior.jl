"""
    MultivariatePoissonProcessPrior{R1,R2}

Gamma prior on all the event rates of a `MultivariatePoissonProcess`.

# Fields

- `λ_α::Vector{R1}`
- `λ_β::R2`
"""
struct MultivariatePoissonProcessPrior{R1<:Real,R2<:Real}
    λ_α::Vector{R1}
    λ_β::R2
end

function DensityInterface.logdensityof(
    prior::MultivariatePoissonProcessPrior, pp::MultivariatePoissonProcess
)
    l = sum(
        logdensityof(Gamma(prior.λα[m], inv(prior.λβ); check_args=false), pp.λ[m]) for
        m in 1:length(pp)
    )
    return l
end
