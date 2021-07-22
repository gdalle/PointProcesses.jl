"""
    MultivariateHawkesProcess

A homegeneous multivariate Hawkes process with integer marks[^Laub_2015].

[^Laub_2015]: Laub, P. J., Taimre, T., and Pollett, P. K. (2015), “Hawkes Processes,” arXiv:1507.02822 [math, q-fin, stat].

# Fields
- `logλ::Vector{Float64}`: logarithms of the base event rates.
- `logα::Matrix{Float64}`: logarithms of the mutual excitation amplitudes.
- `logβ::Matrix{Float64}`: mutual excitation decay rates.

# Examples

```jldoctest
using Random; Random.seed!(63);
logλ = -ones(2)
logα = -ones(2, 2)
logβ = -ones(2, 2)
pp = MultivariateHawkesProcess(logλ, logα, logβ)
h = rand(pp, 0., 1000.)
pp0 = MultivariateHawkesProcess([0., 0., 0.])
pp_est = fit(pp0, h)

# output

MultivariatePoissonProcess([-0.04499736593073609, 0.9850702472689811, 1.9811393903894554])
```
"""
struct MultivariateHawkesProcess <: MultivariatePointProcess
    logλ::Vector{Float64}
    logα::Matrix{Float64}
    logβ::Matrix{Float64}
end

function all_marks(::Type{<:MultivariateHawkesProcess}, θ::Parameter)
    return 1:length(θ.logλ)
end

function intensity(
    ::Type{MultivariateHawkesProcess},
    θ::Parameter,
    h::History{Int},
    t::Float64,
    m::Int,
)
    intens = exp(θ.logλ[m])
    for (tprev, mprev) in zip(h.times, h.marks)
        if tprev >= t
            break
        else
            intens += exp(θ.logα[mprev, m] + exp(-θ.logβ[mprev, m]) * (t - tprev))
        end
    end
    return intens
end

function ground_intensity_bound(
    ::Type{MultivariateHawkesProcess},
    θ::Parameter,
    h::History{Int},
    t::Float64,
)
    λg = ground_intensity(MultivariateHawkesProcess, θ, h, t)
    return λg, λg / 2
end
