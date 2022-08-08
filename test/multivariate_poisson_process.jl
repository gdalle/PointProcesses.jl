using ForwardDiff
using PointProcesses
using Random
using Statistics
using Test
using Zygote

rng = Random.seed!(63)

pp = MultivariatePoissonProcess(rand(rng, 10))

h1 = rand(rng, pp, 0.0, 1000.0)
h2 = simulate_ogata(rng, pp, 0.0, 1000.0)

pp_est1 = fit_mle(MultivariatePoissonProcess{Float32}, h1)
pp_est2 = fit_mle(MultivariatePoissonProcess{Float32}, h2)

λ_error1 = mean(abs, pp_est1.λ - pp.λ)
λ_error2 = mean(abs, pp_est2.λ - pp.λ)

l = logdensityof(pp, h1)
l_est = logdensityof(pp_est1, h1)

f(λ) = logdensityof(MultivariatePoissonProcess(λ), h1)
gf = ForwardDiff.gradient(f, 3 * ones(10))
gz = Zygote.gradient(f, 3 * ones(10))[1]

@test λ_error1 < 0.1
@test λ_error2 < 0.1
@test l_est > l
@test all(gf .≈ gz)
@test all(gf .< 0)
@test all(gz .< 0)
