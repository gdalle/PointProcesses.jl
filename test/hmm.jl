@testset verbose = true "HMM" begin

    @testset "Discrete" begin
        tr = DiscreteMarkovChain([0.3, 0.7], [0.9 0.1; 0.2 0.8])
        em1 = Dists.Normal(1, 0.3)
        em2 = Dists.Normal(2, 0.3)
        hmm = HiddenMarkovModel(tr, [em1, em2])

        states, observations = rand(hmm, 1000)

        tr_init = DiscreteMarkovChain(randprobvec(2), randtransmat(2))
        em1_init = Dists.Normal(rand(), 1)
        em2_init = Dists.Normal(rand(), 1)
        hmm_init = HiddenMarkovModel(tr_init, [em1_init, em2_init])

        hmm_est, logL_evolution = baum_welch(hmm_init, observations, iterations = 100)

        error = transition_matrix(hmm_est) - transition_matrix(hmm)

        @test maximum(error) < 0.1
    end

    @testset "Continuous" begin
        tr = ContinuousMarkovChain([0.3, 0.7], [-1.0 1.0; 2.0 -2.0])
        em1 = PoissonProcess(1.0, Dists.Normal(1, 1))
        em2 = PoissonProcess(2.0, Dists.Normal(-1, 1))
        mmpp = MarkovModulatedPoissonProcess(tr, [em1, em2])

        state_history, observations = rand(mmpp, 0.0, 1000.0)
        nb_events(observations)

    end
end
