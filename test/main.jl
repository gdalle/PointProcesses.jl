using PointProcesses

λ0 = [1.0, 2.0, 0.5]
pp = PoissonProcess(λ0)
history = rand(pp, 0.0, 100.0)
λ_est = learn_poisson(history)