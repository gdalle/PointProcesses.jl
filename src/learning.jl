import Distributions: logpdf, fit
using Flux
using GalacticOptim
using Optim
using Quadrature
using Zygote

ground_intensity_aux(t, p) = ground_intensity(p[1], p[2], t)

function integrated_ground_intensity(pp, history)
    lb, ub, p = get_tmin(history), get_tmax(history), [pp]
    f = (t, p) -> sum(p[1].λ)
    prob = QuadratureProblem(f, lb, ub, p)
    sol = solve(prob, HCubatureJL(), reltol = 1e-3, abstol = 1e-3)
    return sol[1]
end

function integrated_ground_intensity2(λ, history)
    lb, ub, p = get_tmin(history), get_tmax(history), [λ]
    f = (t, p) -> sum(p[1])
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
