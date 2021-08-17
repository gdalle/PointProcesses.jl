# Markov chains

```@meta
DocTestSetup = quote
    using PointProcesses
end
```

```jldoctest markov
julia> using Random

julia> Random.seed!(63);
```

Some point processes are based on underlying Markov processes, which is why we provide a basic implementation for them.

## Discrete time

```jldoctest markov
julia> dmc = DiscreteMarkovChain(π0 = [0.3, 0.7], P = [0.9 0.1; 0.2 0.8]);

julia> states = rand(dmc, 1000);

julia> dmc_est = fit(DiscreteMarkovChain, states);

julia> round.(dmc_est.P, digits=2)
2×2 Matrix{Float64}:
 0.9   0.1
 0.25  0.75
```

### Continuous time

```jldoctest markov
julia> cmc = ContinuousMarkovChain(π0 = [0.3, 0.7], Q = [-1. 1.; 2. -2.]);

julia> history = rand(cmc, 0., 1000.);

julia> cmc_est = fit(ContinuousMarkovChain, history);

julia> round.(cmc_est.Q, digits=2)
2×2 Matrix{Float64}:
 -1.0    1.0
  1.92  -1.92
```