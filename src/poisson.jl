struct PoissonProcess{T} <: PointProcess{Int}
    logλ::Vector{T}
end

PoissonProcess{T}(x::NamedTuple) where {T} = PoissonProcess{T}(x.λ)
PoissonProcess{T}(x::ComponentVector) where {T} = PoissonProcess{T}(x.λ)

function default_params(::Type{PoissonProcess{T}}, history::History{Int}) where {T}
    # return (λ = zeros(T, maximum(history.marks)),)
    return ComponentVector(λ = zeros(T, maximum(history.marks)),)
end

function intensity(pp::PoissonProcess, history::History{Int}, t, m::Int)
    return exp(pp.logλ[m])
end

function mark_distribution(pp::PoissonProcess, history::History{Int}, t)
    return Categorical(exp.(pp.logλ) / sum(exp.(pp.logλ)))
end

function ground_intensity(pp::PoissonProcess, history::History{Int}, t)
    return sum(exp.(pp.logλ))
end

function ground_intensity_bound(pp::PoissonProcess, history::History{Int}, t)
    return ground_intensity(pp, history, t), Inf
end
