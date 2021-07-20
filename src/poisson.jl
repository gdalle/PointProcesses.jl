struct MultivariatePoissonProcess <: PointProcess{Int}
    logλ::Vector{Float64}
end

function default_θ(::Type{MultivariatePoissonProcess}, h)
    return ComponentVector(logλ = zeros(maximum(h.marks)))
end

function intensity(::Type{MultivariatePoissonProcess}, θ, h, t, m)
    return exp(θ.logλ[m])
end

function mark_distribution(::Type{MultivariatePoissonProcess}, θ, h, t)
    return Categorical(exp.(θ.logλ) / sum(exp.(θ.logλ)))
end

function ground_intensity(::Type{MultivariatePoissonProcess}, θ, h, t)
    return sum(exp.(θ.logλ))
end

function ground_intensity_bound(::Type{MultivariatePoissonProcess}, θ, h, t)
    return ground_intensity(MultivariatePoissonProcess, θ, h, t), Inf
end

function integrated_ground_intensity(::Type{MultivariatePoissonProcess}, θ, h)
    return sum(exp.(θ.logλ)) * duration(h)
end
