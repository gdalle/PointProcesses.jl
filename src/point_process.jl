abstract type PointProcess{M} end

function get_θ(pp::PointProcess)
    return ComponentVector{Float64}(ntfromstruct(pp))
end

function intensity(pp::PointProcess{M}, history::History{M}, t, m::M) where {M}
    return intensity(typeof(pp), get_θ(pp), history, t, m)
end

function mark_distribution(pp::PointProcess{M}, history::History{M}, t) where {M}
    return mark_distribution(typeof(pp), get_θ(pp), history, t)
end

function ground_intensity(pp::PointProcess{M}, history::History{M}, t) where {M}
    return ground_intensity(typeof(pp), get_θ(pp), history, t)
end

function ground_intensity_bound(pp::PointProcess{M}, history::History{M}, t) where {M}
    return ground_intensity_bound(typeof(pp), get_θ(pp), history, t)
end

function integrated_ground_intensity(pp::PointProcess{M}, history::History{M}) where {M}
    return integrated_ground_intensity(typeof(pp), get_θ(pp), history)
end

function Distributions.logpdf(pp::PointProcess{M}, history::History{M}) where {M}
    return logpdf(typeof(pp), get_θ(pp), history)
end

## Learning

function integrated_ground_intensity(pptype::Type{<:PointProcess}, θ, history)
    prob = QuadratureProblem(
        (t, θ) -> ground_intensity(pptype, θ, history, t),
        history.tmin,
        history.tmax,
        θ,
    )
    sol = solve(prob, HCubatureJL(), reltol = 1e-3, abstol = 1e-3)
    return sol[1]
end

function Distributions.logpdf(pptype::Type{<:PointProcess}, θ, history)
    l = -integrated_ground_intensity(pptype, θ, history)
    for (t, m) in zip(history.times, history.marks)
        l += log(intensity(pptype, θ, history, t, m))
    end
    return l
end

function Distributions.fit(pptype::Type{<:PointProcess}, history)
    f = OptimizationFunction(
        (θ, params) -> -logpdf(pptype, θ, history),
        GalacticOptim.AutoForwardDiff(),
    )
    prob = OptimizationProblem(f, default_θ(pptype, history))
    sol = solve(prob, LBFGS())
    θ_opt = sol.minimizer
    return θ_opt
end
