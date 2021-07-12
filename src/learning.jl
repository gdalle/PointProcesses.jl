import Distributions: logpdf
using Cuba
using Flux
using GalacticOptim
using Optim
using Quadrature
using Zygote

function logpdf(pp::PoissonProcess, history::History{Int})
    l = 0.
    for (t, m) in zip(get_times(history), get_marks(history))
        l += log(intensity(pp, history, t, m))
    end
    l -= integrated_ground_intensity(pp, history)
    return l
end

function learn_poisson(history::History{Int})
    dim = maximum(get_marks(history))
    f = OptimizationFunction(
        (logλ, p) -> -logpdf(PoissonProcess(exp.(logλ)), history),
        GalacticOptim.AutoZygote()
    )
    prob = OptimizationProblem(f, zeros(dim))
    sol = solve(prob, LBFGS())
    λ = exp.(sol.minimizer)
    return PoissonProcess(λ)
end
