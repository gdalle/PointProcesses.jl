struct PoissonProcess <: PointProcess{Int}
    λ::Float64
    mark_distribution::Categorical
end

function ground_intensity(poissonprocess::PoissonProcess, history, t)
    return poissonprocess.λ
end

function ground_intensity_bound(poissonprocess::PoissonProcess, history, t)
    return poissonprocess.λ
end

function ground_intensity_bound_validity_duration(poissonprocess::PoissonProcess, history, t)
    return Inf
end

function mark_distribution(poissonprocess::PoissonProcess, history, t)
    return poissonprocess.mark_distribution
end