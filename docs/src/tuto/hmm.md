# Hidden Markov models

```@meta
DocTestSetup = quote
    using PointProcesses
end
```

```jldoctest hmm
julia> using Random

julia> using Distributions

julia> Random.seed!(63);
```

## Discrete time

Here is an example of the Baum-Welch estimation algorithm applied to a discrete HMM.

```jldoctest hmm
julia> hmm = HiddenMarkovModel(
           transitions = DiscreteMarkovChain(π0 = [0.3, 0.7], P = [0.9 0.1; 0.2 0.8]),
           emissions = [Normal(1, 0.3), Normal(2, 0.3)]
       );

julia> states, observations = rand(hmm, 1000);

julia> hmm_init = HiddenMarkovModel(
           transitions = DiscreteMarkovChain(π0 = randprobvec(2), P = randtransmat(2)),
           emissions = [Normal(rand(), 1), Normal(rand(), 1)]
       );

julia> hmm_est, logL_evolution = baum_welch(hmm_init, observations, iterations=100);

julia> minimum(diff(logL_evolution)) > -1e-10
true

julia> round.(transition_matrix(hmm_est), digits=2)
2×2 Matrix{Float64}:
 0.9   0.1
 0.23  0.77

julia> round(mean(emission(hmm_est, 1)), digits=2)
1.02

julia> round(mean(emission(hmm_est, 2)), digits=2)
1.99
```

## Continuous time

> Todo