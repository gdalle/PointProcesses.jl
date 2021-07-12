using Distributions

struct PoissonProcess <: PointProcess{Int}
    λ::Vector{Float64}
end

function ground_intensity(pp::PoissonProcess, history::History{Int}, t)
    return sum(pp.λ)
end

function intensity(pp::PoissonProcess, history::History{Int}, t, m::Int)
    return pp.λ[m]
end

function mark_distribution(pp::PoissonProcess, history::History{Int}, t)
    return Categorical(pp.λ / sum(pp.λ))
end

function integrated_ground_intensity(pp::PoissonProcess, history::History{Int})
    return sum(pp.λ) * (get_tmax(history) - get_tmin(history))
end

function ground_intensity_bound(pp::PoissonProcess, history::History{Int}, t)
    return sum(pp.λ)
end

function ground_intensity_bound_validity(pp::PoissonProcess, history::History{Int}, t)
    return Inf
end
