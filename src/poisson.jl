struct PoissonProcess <: PointProcess{Int}
    logλ::Vector{Float64}
end

function default_θ(::Type{PoissonProcess}, history)
    return ComponentVector(λ = zeros(maximum(history.marks)))
end

function intensity(::Type{PoissonProcess}, θ, history, t, m)
    return exp(θ.logλ[m])
end

function mark_distribution(::Type{PoissonProcess}, θ, history, t)
    return Categorical(exp.(θ.logλ) / sum(exp.(θ.logλ)))
end

function ground_intensity(::Type{PoissonProcess}, θ, history, t)
    return sum(exp.(θ.logλ))
end

function ground_intensity_bound(::Type{PoissonProcess}, θ, history, t)
    return ground_intensity(PoissonProcess, θ, history, t), Inf
end
