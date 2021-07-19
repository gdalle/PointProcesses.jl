using PointProcesses
using GalacticOptim
using Optim
using Zygote
using Quadrature
using ForwardDiff
using ComponentArrays

logλ0 = [-2, 0, 2]
pp = PoissonProcess(logλ0)
θ0 = get_θ(pp)

history = rand(pp, 0.0, 10.0)

integrated_ground_intensity(PoissonProcess, θ0, history)
logpdf(PoissonProcess, θ0, history)

g = ForwardDiff.gradient(θ -> integrated_ground_intensity(PoissonProcess, θ, history), θ0)
g = Zygote.gradient(θ -> integrated_ground_intensity(PoissonProcess, θ, history), θ0)

θ_init = ComponentVector(logλ = [0., 0., 0.])
fit(PoissonProcess, θ_init, history)
