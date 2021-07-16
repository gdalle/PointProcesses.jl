using PointProcesses
using GalacticOptim
using Optim
using Zygote
using Quadrature
using ForwardDiff

λ0 = [1.0, 2.0, 0.5]
pp = PoissonProcess(λ0)
history = rand(pp, 0.0, 10.0)

integrated_ground_intensity(PoissonProcess(λ0), history)
logpdf(pp, history)

g = ForwardDiff.gradient(λ -> integrated_ground_intensity(PoissonProcess(λ), history), ones(3))