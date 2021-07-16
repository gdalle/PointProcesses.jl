using PointProcesses
using GalacticOptim
using Optim
using Zygote
using Quadrature

λ0 = [1.0, 2.0, 0.5]
pp = PoissonProcess(λ0)
history = rand(pp, 0.0, 10.0)

integrated_ground_intensity2(λ0, history)
integrated_ground_intensity(PoissonProcess(λ0), history)

logpdf(pp, history)

Zygote.gradient(λ -> integrated_ground_intensity2(λ, history), ones(3))
Zygote.gradient(λ -> integrated_ground_intensity(PoissonProcess(λ), history), ones(3))

fit(PoissonProcess, history)