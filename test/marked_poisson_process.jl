pp = MarkedPoissonProcess(1.0, Categorical([0.1, 0.3, 0.6]))
h = rand(pp, 0.0, 1000.0)
pp_est = fit_mle(MarkedPoissonProcess{Float32,Int,Categorical}, h)
λ_error = mean(abs, pp_est.λ - pp.λ)

@test λ_error < 0.1
