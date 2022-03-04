@testset verbose = true "Models" begin
    @testset "Multivariate Poisson" begin
        pp = MultivariatePoissonProcess(rand(10))
        h = rand(pp, 0.0, 1000.0)
        pp_est = fit_mle(MultivariatePoissonProcess, h)
        error = pp_est.λ - pp.λ
        @test maximum(error) < 0.2
    end

    @testset "Generic Poisson" begin
    end
end
