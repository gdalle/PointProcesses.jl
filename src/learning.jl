import Distributions: logpdf, fit
using ForwardDiff
using GalacticOptim
using Optim
using Quadrature

ground_intensity_aux(t, params) = ground_intensity(params[1], params[2], t)

function integrated_ground_intensity(pp, history)
    lb, ub = get_tmin(history), get_tmax(history)
    params = [pp, history]
    prob = QuadratureProblem(ground_intensity_aux, lb, ub, params)
    sol = solve(prob, HCubatureJL(), reltol = 1e-3, abstol = 1e-3)
    return sol[1]
end

function logpdf(pp::PointProcess, history)
    l = 0.0
    for (t, m) in zip(get_times(history), get_marks(history))
        l += log(intensity(pp, history, t, m))
    end
    l -= integrated_ground_intensity(pp, history)
    return l
end

neglogpdf_aux(x, params) = -logpdf(params[1](x), params[2])

function fit(pptype::Type{<:PointProcess}, history)
    x0 = default_param(pptype, history)
    params = [pptype, history]
    prob = OptimizationProblem(f, x0, params)
    sol = solve(prob, LBFGS())
    x_opt = sol.minimizer
    return pptype(x_opt)
end
