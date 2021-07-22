"""
    MultivariatePoissonProcess

A homegeneous Poisson process with contiguous integer marks.

# Fields
- `logλ::Vector{Float64}`: logarithms of the event rates.

# Examples

```jldoctest
using Random; Random.seed!(63);
pp = MultivariatePoissonProcess([0., 1., 2.])
h = rand(pp, 0., 1000.)
pp0 = MultivariatePoissonProcess([0., 0., 0.])
pp_est = fit(pp0, h)

# output

MultivariatePoissonProcess([-0.04499736593073609, 0.9850702472689811, 1.9811393903894554])
```
"""
struct MultivariatePoissonProcess <: MultivariatePointProcess
    logλ::Vector{Float64}
end

function all_marks(::Type{<:MultivariatePoissonProcess}, θ::Parameter)
    return 1:length(θ.logλ)
end

function intensity(
    ::Type{MultivariatePoissonProcess},
    θ::Parameter,
    h::History{Int},
    t::Float64,
    m::Int,
)
    return exp(θ.logλ[m])
end

function ground_intensity_bound(
    ::Type{MultivariatePoissonProcess},
    θ::Parameter,
    h::History{Int},
    t::Float64,
)
    λg = ground_intensity(MultivariatePoissonProcess, θ, h, t)
    return λg, Inf
end
