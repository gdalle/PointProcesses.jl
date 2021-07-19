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

    history = rand(pp, 0.0, 10000.0)

    integrated_ground_intensity(pp, history)

    logpdf(pp, history)

    g = ForwardDiff.gradient(θ -> logpdf(MultivariatePoissonProcess, θ, history), θ0)
    g = Zygote.gradient(θ -> logpdf(MultivariatePoissonProcess, θ, history), θ0)

    θ_est = fit(MultivariatePoissonProcess, history)
    @test maximum(abs.(θ_est.logλ - θ0.logλ)) < 1e-1
end
