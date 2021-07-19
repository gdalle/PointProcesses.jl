using PointProcesses
using GalacticOptim
using Optim
using Zygote
using Quadrature
using ForwardDiff
using ParameterHandling

λ0 = [-2, 0, 2]
pp = PoissonProcess(λ0)
history = rand(pp, 0.0, 10.0)

integrated_ground_intensity(PoissonProcess(λ0), history)
logpdf(pp, history)

g = ForwardDiff.gradient(
    λ -> integrated_ground_intensity(PoissonProcess(λ), history),
    ones(3),
)

g = Zygote.gradient(
    λ -> integrated_ground_intensity(PoissonProcess(λ), history),
    ones(3),
)

pptype = PoissonProcess{Float64}

x0 = default_params(pptype, history)
x0, unpack = flatten(x0_tuple)
f = OptimizationFunction(
    (x, params) -> -logpdf(pptype(unpack(x)), history),
    GalacticOptim.AutoForwardDiff(),
)
prob = OptimizationProblem(f, x0)
sol = solve(prob, LBFGS())
x_opt = sol.minimizer
pp_opt = PoissonProcess{Float64}(unpack(x_opt))

fit(PoissonProcess{Float64}, history)
