ground_intensity_aux(t, θ) = ground_intensity(θ[1], θ[2], t)

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

function Distributions.fit(pptype::Type{<:PointProcess}, initial_θ, history)
    f = OptimizationFunction(
        (θ, params) -> -logpdf(pptype, θ, history),
        GalacticOptim.AutoZygote(),
    )
    prob = OptimizationProblem(f, initial_θ)
    sol = solve(prob, LBFGS())
    θ_opt = sol.minimizer
    return θ_opt
end
