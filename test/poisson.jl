@testset "Poisson processes" begin
    M = 100
    λ = rand(M)
    λ_init = 0.5 * ones(M)

    pp1 = MultivariatePoissonProcess(λ)
    pp2 = NaiveMultivariatePoissonProcess(λ)
    pp2_init = NaiveMultivariatePoissonProcess(λ_init)

    history1 = rand(pp1, 0., 100.)
    history2 = rand(pp2, 0., 100.)

    pp1_est = fit(MultivariatePoissonProcess, history1)
    pp2_est = fit(pp2_init, history2)

    @test sum(abs.(pp1_est.λ .- pp1.λ)) / M < 1e-1
    @test sum(abs.(pp2_est.λ .- pp2.λ)) / M < 1e-1
end
