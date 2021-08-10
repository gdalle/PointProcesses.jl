# Tutorial

```@meta
DocTestSetup = quote
    using PointProcesses
end
```

## Working with event histories

```jldoctest
julia> h = TemporalHistory([0.2, 0.8, 1.1], ["a", "b", "c"], 0.0, 2.0)
TemporalHistory{String}([0.2, 0.8, 1.1], ["a", "b", "c"], 0.0, 2.0)

julia> duration(h)
2.0

julia> nb_events(h)
3

julia> nb_events(h, 1.0, 2.0)
1

julia> has_events(h)
true

julia> has_events(h, 1.5, 2.0)
false

julia> push!(h, 1.7, "d")

julia> has_events(h, 1.5, 2.0)
true
```

## Working with Markov chains

### Discrete time

```jldoctest
using Random; Random.seed!(63)
dmc = DiscreteMarkovChain([0.3, 0.7], [0.9 0.1; 0.2 0.8])
states = rand(dmc, 100)
dmc_est = fit(DiscreteMarkovChain, states)
round.(dmc_est.P, digits=3)

# output

2×2 Matrix{Float64}:
 0.943  0.057
 0.417  0.583
```

### Continuous time

```jldoctest
using Random; Random.seed!(63)
cmc = ContinuousMarkovChain([0.3, 0.7], [-1. 1.; 2. -2.])
h = rand(cmc, 0., 100.)
cmc_est = fit(ContinuousMarkovChain, h)
round.(cmc_est.Q, digits=3)

# output

2×2 Matrix{Float64}:
 -1.175   1.175
  1.793  -1.793
```

## Working with Poisson processes

```jldoctest
using Random; Random.seed!(63);
pp = TemporalPoissonProcess([0.5, 1., 2.])
h = rand(pp, 0., 100.)
pp_init = TemporalPoissonProcess([1., 1., 1.])
pp_est = fit(pp_init, h)
round.(pp_est.λ, digits=3)

# output

3-element Vector{Float64}:
 0.6
 1.14
 1.79
```

## Working with Hidden Markov Models

```jldoctest
using Random; Random.seed!(63)
dmc = DiscreteMarkovChain([0.3, 0.7], [0.9 0.1; 0.2 0.8])
emission1 = BoundedTemporalPointProcess(TemporalPoissonProcess([0., 1., 2.]), 0., 1.)
emission2 = BoundedTemporalPointProcess(TemporalPoissonProcess([2., 1., 0.]), 0., 1.)
hmm = HiddenMarkovModel(dmc, [emission1, emission2])
states, observations = rand(hmm, 100)
sum(nb_events(observations[t]) for t = 1:100)

# output

268
```