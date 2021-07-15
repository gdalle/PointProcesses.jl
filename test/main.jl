using PointProcesses
using GalacticOptim
using Optim
using Zygote
using Quadrature

λ0 = [1.0, 2.0, 0.5]
pp = PoissonProcess(λ0)
history = rand(pp, 0.0, 10.0)

integrated_ground_intensity(pp, history)
logpdf(pp, history)

pptype = PoissonProcess
x0 = default_param(pptype, history)
Zygote.gradient(x -> integrated_ground_intensity(pptype(x), history), ones(3))

fit(PoissonProcess, history)