@testset "Poisson processes" begin
    M = 100
    λ = rand(M)
    λ_init = 0.5 * ones(M)

    pp1 = PoissonProcess(λ)
    pp2 = NaivePoissonProcess(λ)
    pp2_init = NaivePoissonProcess(λ_init)

    history1 = rand(pp1, 0., 100.)
    history2 = rand(pp2, 0., 100.)

    pp1_est = fit(PoissonProcess, history1)
    pp2_est = fit(pp2_init, history2)

    @test sum(abs.(pp1_est.λ .- pp1.λ)) / M < 1e-1
    @test sum(abs.(pp2_est.λ .- pp2.λ)) / M < 1e-1
end
