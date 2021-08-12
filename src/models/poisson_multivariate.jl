"""
    MultivariatePoissonProcess{R}

Homogeneous temporal multivariate Poisson process.

# Fields
- `λ::Vector{R}`: event rates.
"""
@with_kw struct MultivariatePoissonProcess{R} <: TemporalPointProcess{Int}
    λ::Vector{R}
end

function Base.rand(rng::AbstractRNG, pp::MultivariatePoissonProcess, tmin, tmax)
    N = rand(Poisson(sum(pp.λ) * (tmax - tmin)))
    times = rand(Uniform(tmin, tmax), N)
    marks = rand(Categorical(pp.λ / sum(pp.λ)), N)
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
