@testset verbose = true "Markov" begin

    @testset "Discrete" begin
        dmc = DiscreteMarkovChain([0.3, 0.7], [0.9 0.1; 0.2 0.8])
        states = rand(dmc, 1000)
        dmc_est = fit(DiscreteMarkovChain, states)
        error = transition_matrix(dmc_est) - transition_matrix(dmc)
        @test maximum(error) < 0.2
    end

    @testset "Continuous" begin
        cmc = ContinuousMarkovChain([0.3, 0.7], [-1. 1.; 2. -2.])
        history = rand(cmc, 0., 1000.)
        cmc_est = fit(ContinuousMarkovChain, history)
        error = rate_matrix(cmc_est) - rate_matrix(cmc)
        @test maximum(error) < 0.2
    end

end
