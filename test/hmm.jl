@testset verbose = true "HMM" begin
    tr = DiscreteMarkovChain([0.3, 0.7], [0.9 0.1; 0.2 0.8])
    em1 = Normal(1, 0.3)
    em2 = Normal(2, 0.6)
    hmm = HiddenMarkovModel(tr, [em1, em2])

    states, observations = rand(hmm, 1000)

    tr_init = DiscreteMarkovChain(randprobvec(2), randtransmat(2))
    em1_init = Normal(rand(), 1)
    em2_init = Normal(rand(), 1)
    hmm_init = HiddenMarkovModel(tr_init, [em1_init, em2_init])

    hmm_est1, logL_evolution1 = baum_welch(hmm_init, observations, iterations = 100)
    hmm_est2, logL_evolution2 = baum_welch(hmm_init, observations, iterations = 100, log=true)

    transition_error1 = abs.(transition_matrix(hmm_est1) - transition_matrix(hmm))
    transition_error2 = abs.(transition_matrix(hmm_est2) - transition_matrix(hmm))

    @test maximum(transition_error1) < 0.1
    @test maximum(transition_error2) < 0.1

    μ_error1 = abs.([emission(hmm_est1, s).μ - emission(hmm, s).μ for s = 1:2])
    μ_error2 = abs.([emission(hmm_est2, s).μ - emission(hmm, s).μ for s = 1:2])

    @test maximum(μ_error1) < 0.1
    @test maximum(μ_error2) < 0.1

    σ_error1 = abs.([emission(hmm_est1, s).σ - emission(hmm, s).σ for s = 1:2])
    σ_error2 = abs.([emission(hmm_est2, s).σ - emission(hmm, s).σ for s = 1:2])

    @test maximum(σ_error1) < 0.1
    @test maximum(σ_error2) < 0.1
end
