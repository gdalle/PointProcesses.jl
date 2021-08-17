"""
    MultivariatePoissonProcess{R}

Homogeneous temporal multivariate Poisson process.

# Fields
- `λ::Vector{R}`: event rates.
"""
struct MultivariatePoissonProcess{R} <: TemporalPointProcess{Int}
    λ::Vector{R}
end

function Base.rand(rng::AbstractRNG, pp::MultivariatePoissonProcess, tmin, tmax)
    N = rand(rng, Poisson(sum(pp.λ) * (tmax - tmin)))
    times = rand(rng, Uniform(tmin, tmax), N)
    marks = rand(rng, Categorical(pp.λ / sum(pp.λ)), N)
    return TemporalHistory(times, marks, tmin, tmax)
end

function Distributions.fit(::Type{MultivariatePoissonProcess}, history::TemporalHistory)
    M = maximum(history.marks)
    λ = zeros(Float64, M)
    for m in history.marks
        λ[m] += 1.
    end
    λ ./= duration(history)
    return MultivariatePoissonProcess(λ)
end
