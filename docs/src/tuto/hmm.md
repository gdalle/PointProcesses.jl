# Hidden Markov models

```@meta
DocTestSetup = quote
    using PointProcesses
end
```

## Discrete time

Here is an example of the Baum-Welch estimation algorithm applied to a discrete HMM.

```jldoctest hmm
using Random, MeasureTheory; Random.seed!(63)

tr = DiscreteMarkovChain([0.3, 0.7], [0.9 0.1; 0.2 0.8])
em1 = Dists.Normal(1, 0.3)
em2 = Dists.Normal(2, 0.3)
hmm = HiddenMarkovModel(tr, [em1, em2])

states, observations = rand(hmm, 1000)

tr_init = DiscreteMarkovChain(randprobvec(2), randtransmat(2))
em1_init = Dists.Normal(rand(), 1)
em2_init = Dists.Normal(rand(), 1)
hmm_init = HiddenMarkovModel(tr_init, [em1_init, em2_init])

hmm_est, logL_evolution = baum_welch(hmm_init, observations, iterations=100)

error = transition_matrix(hmm_est) - transition_matrix(hmm)

maximum(error) < 0.1

# output

true
```

## Continuous time

Here is an example of the Baum-Welch estimation algorithm applied to a MMPP.

```jldoctest mmpp
using Random, MeasureTheory; Random.seed!(63)

tr = ContinuousMarkovChain([0.3, 0.7], [-1. 1.; 2. -2.])
em1 = PoissonProcess(1., Dists.Normal(1, 1))
em2 = PoissonProcess(2., Dists.Normal(-1, 1))
mmpp = MarkovModulatedPoissonProcess(tr, [em1, em2])

state_history, observations = rand(mmpp, 0., 1000.)
nb_events(observations)

# output

1345
```