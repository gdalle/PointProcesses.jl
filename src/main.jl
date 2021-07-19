using ForwardDiff
using GalacticOptim
using Optim
using PointProcesses
using Quadrature
using Zygote

pp = PoissonProcess([-2, 0, 2])
θ0 = get_θ(pp)

history = rand(pp, 0.0, 1000.0)

integrated_ground_intensity(pp, history)

logpdf(pp, history)

g = ForwardDiff.gradient(θ -> logpdf(PoissonProcess, θ, history), θ0)
g = Zygote.gradient(θ -> logpdf(PoissonProcess, θ, history), θ0)

@time fit(PoissonProcess, history)
