using Documenter
using ForwardDiff
using GalacticOptim
using Optim
using PointProcesses
using Quadrature
using Test
using Zygote

@testset "PointProcesses.jl" begin
    pp = MultivariatePoissonProcess([-2, 0, 2])
    θ0 = get_θ(pp)
    h = rand(pp, 0.0, 10000.0)
    integrated_ground_intensity(pp, h)
    logpdf(pp, h)
    g = ForwardDiff.gradient(θ -> logpdf(MultivariatePoissonProcess, θ, h), θ0)
    θ_est = fit(MultivariatePoissonProcess, h)
    @test maximum(abs.(θ_est.logλ - θ0.logλ)) < 1e-1
end
