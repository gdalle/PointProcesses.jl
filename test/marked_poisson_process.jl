using DensityInterface
using Distributions
using ForwardDiff
using PointProcesses
using Statistics
using StatsAPI
using Test
# using Zygote

rng = Random.seed!(63)

pp = MarkedPoissonProcess(1.0, Categorical([0.1, 0.3, 0.6]))

h1 = rand(rng, pp, 0.0, 1000.0)
h2 = simulate_ogata(rng, pp, 0.0, 1000.0)

pp_est1 = fit(MarkedPoissonProcess{Int,Float32,Categorical}, [h1, h1])
pp_est2 = fit(MarkedPoissonProcess{Int,Float32,Categorical}, [h2, h2])

λ_error1 = mean(abs, pp_est1.λ - pp.λ)
λ_error2 = mean(abs, pp_est2.λ - pp.λ)
p_error1 = mean(abs, pp_est1.mark_dist.p - pp.mark_dist.p)
p_error2 = mean(abs, pp_est2.mark_dist.p - pp.mark_dist.p)

l = logdensityof(pp, h1)
l_est = logdensityof(pp_est1, h1)

f2(λ) = logdensityof(MarkedPoissonProcess(λ, Categorical([0.1, 0.3, 0.6])), h1)
gf = ForwardDiff.derivative(f, 3)
# gz = Zygote.gradient(f, 3)[1]

@test DensityKind(pp) == HasDensity()
@test λ_error1 < 0.1
@test λ_error2 < 0.1
@test p_error1 < 0.1
@test p_error2 < 0.1
@test l_est > l
# @test all(gf .≈ gz)
@test all(gf .< 0)
# @test all(gz .< 0)
