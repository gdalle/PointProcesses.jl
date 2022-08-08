pp = MultivariatePoissonProcess(rand(10))
h = rand(pp, 0.0, 1000.0)
pp_est = fit_mle(MultivariatePoissonProcess{Float32}, h)
error = mean(abs, pp_est.λ - pp.λ)

@test error < 0.1
