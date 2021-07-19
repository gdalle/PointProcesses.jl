ground_intensity_aux(t, params) = ground_intensity(params[1], params[2], t)

function integrated_ground_intensity(pp, history)
    params = [pp, history]
    prob = QuadratureProblem(ground_intensity_aux, history.tmin, history.tmax, params)
    sol = solve(prob, HCubatureJL(), reltol = 1e-3, abstol = 1e-3)
    return sol[1]
end

function Distributions.logpdf(pp::PointProcess, history)
    l = -integrated_ground_intensity(pp, history)
    for (t, m) in zip(history.times, history.marks)
        l += log(intensity(pp, history, t, m))
    end
    return l
end

function Distributions.fit(pptype::Type{<:PointProcess}, history)
    x0_tuple = default_params(pptype, history)
    x0, unpack = flatten(x0_tuple)
    println(x0_tuple, " ", x0)
    f = OptimizationFunction(
        (x, params) -> -logpdf(pptype(unpack(x)), history),
        GalacticOptim.AutoZygote(),
    )
    prob = OptimizationProblem(f, x0)
    sol = solve(prob, LBFGS())
    x_opt = sol.minimizer
    return pptype(unpack(x_opt))
end
