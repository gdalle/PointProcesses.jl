# Markov chains

```@meta
DocTestSetup = quote
    using PointProcesses
end
```

Some point processes are based on underlying Markov processes, which is why we provide a basic implementation for them.

## Discrete time

```jldoctest
using Random; Random.seed!(63)

dmc = DiscreteMarkovChain([0.3, 0.7], [0.9 0.1; 0.2 0.8])
states = rand(dmc, 1000)
dmc_est = fit(DiscreteMarkovChain, states)
round.(dmc_est.P, digits=2)

# output

2×2 Matrix{Float64}:
 0.9   0.1
 0.25  0.75
```

## Continuous time

```jldoctest
using Random; Random.seed!(63)

cmc = ContinuousMarkovChain([0.3, 0.7], [-1. 1.; 2. -2.])
history = rand(cmc, 0., 1000.)
cmc_est = fit(ContinuousMarkovChain, history)
round.(cmc_est.Q, digits=2)

# output

2×2 Matrix{Float64}:
 -0.93   0.93
  1.87  -1.87
```