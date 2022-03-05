@testset "MMPP" begin
    tr = ContinuousMarkovChain([0.3, 0.7], [-1.0 1.0; 2.0 -2.0])
    em1 = PoissonProcess(1.0, Dists.Normal(1, 1))
    em2 = PoissonProcess(2.0, Dists.Normal(-1, 1))
    mmpp = MarkovModulatedPoissonProcess(tr, [em1, em2])

    state_history, observations = rand(mmpp, 0.0, 1000.0)
    nb_events(observations)
end
