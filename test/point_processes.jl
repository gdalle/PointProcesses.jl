@testset verbose = true "Point processes" begin
    @testset "Interface" begin
        pp = NaiveMultivariatePoissonProcess([1., 2., 3.])
        h = rand(pp, 0., 100.)
        pp_init = NaiveMultivariatePoissonProcess(ones(3))
        pp_est1 = fit(pp_init, h)
        error1 = intensity(pp_est1) - intensity(pp)
        @test maximum(error1) < 0.1
        pp_est2 = fit(pp_init, h, adtype=GalacticOptim.AutoForwardDiff(), alg=Optim.LBFGS())
        error2 = intensity(pp_est2) - intensity(pp)
        @test maximum(error2) < 0.1
    end
end
