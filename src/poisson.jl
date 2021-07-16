using ComponentArrays
using Distributions

struct PoissonProcess{T} <: PointProcess{Int}
    λ::Vector{T}
    PoissonProcess(logλ::Vector{T}) where {T} = new{T}(exp.(logλ))
end

function intensity(pp::PoissonProcess, history::History{Int}, t, m::Int)
    return pp.λ[m]
end

function all_marks(pp::PoissonProcess)
    return 1:length(pp.logλ)
end

function mark_distribution(pp::PoissonProcess, history::History{Int}, t)
    return Categorical(pp.λ / sum(pp.λ))
end

function ground_intensity_bound(pp::PoissonProcess, history::History{Int}, t)
    return ground_intensity(pp, history, t), Inf
end

function default_param(::Type{PoissonProcess{T}}, history::History{Int}) where {T}
    return zeros(T, maximum(get_marks(history)))
end

# function integrated_ground_intensity(pp::PoissonProcess, history::History{Int})
#     return sum(pp.λ) * (get_tmax(history) - get_tmin(history))
# end