@testset verbose = true "Models" begin
    @testset "Poisson" begin
        mark_dist = Dists.Categorical([0.1, 0.3, 0.6])
        pp = PoissonProcess(5.0, mark_dist)
        h = rand(pp, 0.0, 100.0)
        pp_est = fit(PoissonProcess{Dists.Categorical}, h)
        error = Dists.probs(mark_distribution(pp_est)) - Dists.probs(mark_distribution(pp))
        @test maximum(error) < 0.1
    end
end
