using DensityInterface
using ForwardDiff
using PointProcesses
using Random
using Statistics
using StatsAPI
using Test
# using Zygote

rng = Random.seed!(63)

λ = rand(rng, 10)
pp_mark = PoissonProcess(λ)
pp_multi = MultivariatePoissonProcess(λ)
@test ground_intensity(pp_mark) == ground_intensity(pp_multi)
@test mark_distribution(pp_mark) == mark_distribution(pp_multi)

λ = rand(rng, 10)
bpp_mark = BoundedPointProcess(PoissonProcess(λ), 0.0, 1000.0)
bpp_multi = BoundedPointProcess(MultivariatePoissonProcess(λ), 0.0, 1000.0)
@test ground_intensity(bpp_mark, 0, []) == ground_intensity(bpp_multi, 0, [])
@test mark_distribution(bpp_mark, 0, []) == mark_distribution(bpp_multi, 0, [])

h1_mark = rand(Random.seed!(63), pp_mark, 0.0, 1000.0)
h1_multi = rand(Random.seed!(63), pp_multi, 0.0, 1000.0)
@test event_times(h1_mark) == event_times(h1_multi)
h2_mark = simulate_ogata(Random.seed!(63), pp_mark, 0.0, 1000.0)
h2_multi = simulate_ogata(Random.seed!(63), pp_multi, 0.0, 1000.0)
@test event_times(h2_mark) == event_times(h2_multi)
h2bis_mark = rand(Random.seed!(63), bpp_mark)
h2bis_multi = rand(Random.seed!(63), bpp_multi)
@test event_times(h2bis_mark) == event_times(h2bis_multi)

pp_est1_mark = fit(
    PoissonProcess{Float32,Categorical{Float32,Vector{Float32}}}, [h1_mark, h1_mark]
)
pp_est1_multi = fit(MultivariatePoissonProcess{Float32}, [h1_multi, h1_multi])
@test ground_intensity(pp_est1_mark) ≈ ground_intensity(pp_est1_multi)
@test mark_distribution(pp_est1_mark) ≈ mark_distribution(pp_est1_multi)

pp_est2_mark = fit(PoissonProcess{Float32,Categorical}, [h2_mark, h2_mark])
pp_est2_multi = fit(MultivariatePoissonProcess{Float32}, [h2_multi, h2_multi])
@test ground_intensity(pp_est2_mark) ≈ ground_intensity(pp_est2_multi)
@test mark_distribution(pp_est2_mark) ≈ mark_distribution(pp_est2_multi)

α = rand(Uniform(0.5, 1.5), 10)
β = rand()
prior_mark = PoissonProcessPrior(α, β)
pp_est3_mark = fit_map(
    PoissonProcess{Float32,Categorical{Float32,Vector{Float32}}},
    prior_mark,
    [h1_mark, h2_mark],
)
prior_multi = MultivariatePoissonProcessPrior(α, β)
pp_est3_multi = fit_map(
    MultivariatePoissonProcess{Float32}, prior_multi, [h1_multi, h2_multi]
)
@test ground_intensity(pp_est3_mark) ≈ ground_intensity(pp_est3_multi)
@test mark_distribution(pp_est3_mark) ≈ mark_distribution(pp_est3_multi)

@test logdensityof(pp_mark, h1_mark) == logdensityof(pp_multi, h1_mark)

f_mark(λ) = logdensityof(PoissonProcess(λ), h1_mark)
f_multi(λ) = logdensityof(MultivariatePoissonProcess(λ), h1_multi)
@test ForwardDiff.gradient(f_mark, 3 * ones(10)) ==
    ForwardDiff.gradient(f_multi, 3 * ones(10))
