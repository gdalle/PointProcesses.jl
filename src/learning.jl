using Cuba
import Distributions: logpdf, fit
using Flux
using GalacticOptim
using Optim
using Quadrature
using Zygote

function integrated_ground_intensity(pp, history)
    f = (t, p) -> ground_intensity(p[1], p[2], t)
    lb, ub, p = get_tmin(history), get_tmax(history), [pp, history]
    prob = QuadratureProblem(f, lb, ub, p)
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

function fit(pptype::Type{<:PointProcess}, history)
    f = OptimizationFunction((x, p) -> -logpdf(pptype(x), p[1]), GalacticOptim.AutoZygote())
    x0 = default_param(pptype, history)
    p = [history]
    prob = OptimizationProblem(f, x0, p)
    sol = solve(prob, LBFGS())
    x_opt = sol.minimizer
    return pptype(x_opt)
end
