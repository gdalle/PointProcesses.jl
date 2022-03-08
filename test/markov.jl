@testset verbose = true "Markov" begin

    @testset "Discrete" begin
        dmc = DiscreteMarkovChain(; π0=[0.3, 0.7], P=[0.9 0.1; 0.2 0.8])
        @test nb_states(dmc) == 2
        @test stationary_distribution(dmc) ≈ [0.2 / (0.1 + 0.2), 0.1 / (0.1 + 0.2)]

        states = rand(dmc, 1000)
        dmc_est_mle = fit_mle(DiscreteMarkovChain, states)
        error_mle = mean(abs, transition_matrix(dmc_est_mle) - transition_matrix(dmc))
        @test error_mle < 0.1

        prior = DiscreteMarkovChainPrior(; π0α=[1.0, 1.0], Pα=1000 * ones(2, 2))
        dmc_est_map = fit_map(DiscreteMarkovChain, prior, states)
        error_map = mean(abs, transition_matrix(dmc_est_map) - transition_matrix(dmc))
        @test error_map > error_mle
    end

    @testset "Continuous" begin
        cmc = ContinuousMarkovChain(; π0=[0.3, 0.7], Q=[-1.0 1.0; 2.0 -2.0])
        history = rand(cmc, 0.0, 1000.0)
        cmc_est = fit_mle(ContinuousMarkovChain, history)
        error = mean(abs, rate_matrix(cmc_est) - rate_matrix(cmc))
        @test error < 0.1
    end

end
