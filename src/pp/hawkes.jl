"""
    MultivariateHawkesProcess

A homegeneous multivariate Hawkes process with integer marks[^Laub_2015].

# Fields
- `logλ::Vector{Float64}`: logarithms of the base event rates.
- `logα::Matrix{Float64}`: logarithms of the mutual excitation amplitudes.
- `β::Matrix{Float64}`: mutual excitation decay rates.

[^Laub_2015]: Laub, P. J., Taimre, T., and Pollett, P. K. (2015), “Hawkes Processes,” arXiv:1507.02822 [math, q-fin, stat].
"""
struct MultivariateHawkesProcess <: PointProcess{Int}
    logλ::Vector{Float64}
    logα::Matrix{Float64}
    β::Matrix{Float64}
end
