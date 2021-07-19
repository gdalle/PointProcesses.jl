using ForwardDiff
using GalacticOptim
using Optim
using PointProcesses
using Quadrature
using Zygote

pp = MultivariatePoissonProcess([-2, 0, 2])
θ0 = get_θ(pp)

history = rand(pp, 0.0, 10000.0)

integrated_ground_intensity(pp, history)

logpdf(pp, history)

g = ForwardDiff.gradient(θ -> logpdf(MultivariatePoissonProcess, θ, history), θ0)
g = Zygote.gradient(θ -> logpdf(MultivariatePoissonProcess, θ, history), θ0)

@time fit(MultivariatePoissonProcess, history)
