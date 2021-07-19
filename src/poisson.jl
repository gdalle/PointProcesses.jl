struct MultivariatePoissonProcess <: PointProcess{Int}
    logλ::Vector{Float64}
end

function default_θ(::Type{MultivariatePoissonProcess}, history)
    return ComponentVector(logλ = zeros(maximum(history.marks)))
end

function intensity(::Type{MultivariatePoissonProcess}, θ, history, t, m)
    return exp(θ.logλ[m])
end

function mark_distribution(::Type{MultivariatePoissonProcess}, θ, history, t)
    return Categorical(exp.(θ.logλ) / sum(exp.(θ.logλ)))
end

function ground_intensity(::Type{MultivariatePoissonProcess}, θ, history, t)
    return sum(exp.(θ.logλ))
end

function ground_intensity_bound(::Type{MultivariatePoissonProcess}, θ, history, t)
    return ground_intensity(MultivariatePoissonProcess, θ, history, t), Inf
end

function integrated_ground_intensity(::Type{MultivariatePoissonProcess}, θ, history)
    return sum(exp.(θ.logλ)) * duration(history)
end
