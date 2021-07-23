"""
    PoissonProcess{R}

A homogeneous temporal multivariate Poisson process.

# Fields
- `λ::Vector{R}`: event rates.

# Examples

```jldoctest
using Random; Random.seed!(63);
pp = PoissonProcess([0.5, 1., 2.])
h = rand(pp, 0., 100.)
pp_init = PoissonProcess([1., 1., 1.])
pp_est = fit(pp_init, h)

# output

PoissonProcess{Float64}([0.5999999999996618, 1.1400000000005681, 1.7900000000002536])
```
"""
struct PoissonProcess{R} <: MultivariateTemporalPointProcess
    λ::Vector{R}
end

PoissonProcess{R}(nt::NamedTuple) where {R} = PoissonProcess(nt.λ)

function build_transform(pp::PoissonProcess)
    M = length(pp.λ)
    return as((λ = as(Vector, asℝ₊, M),))
end

function all_marks(pp::PoissonProcess)
    return 1:length(pp.λ)
end

function intensity(pp::PoissonProcess, h::History{Int}, t::Float64, m::Int)
    return pp.λ[m]
end

function ground_intensity_bound(pp::PoissonProcess, h::History{Int}, t::Float64)
    λg = ground_intensity(pp, h, t)
    return λg, Inf
end
