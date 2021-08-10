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
fit(DiscreteMarkovChain, states)

# output

DiscreteMarkovChain([1.0, 0.0], [0.9425287356321839 0.05747126436781609; 0.4166666666666667 0.5833333333333334])
```

### Continuous time

```jldoctest
using Random; Random.seed!(63)
cmc = ContinuousMarkovChain([0.3, 0.7], [-1. 1.; 2. -2.])
h = rand(cmc, 0., 100.)
fit(ContinuousMarkovChain, h)

# output

ContinuousMarkovChain([1.0, 0.0], [-1.1753044442034237 1.1753044442034237; 1.7933767144927923 -1.7933767144927923])
```

## Working with Poisson processes

```jldoctest
using Random; Random.seed!(63);
pp = TemporalPoissonProcess([0.5, 1., 2.])
h = rand(pp, 0., 100.)
pp_init = TemporalPoissonProcess([1., 1., 1.])
pp_est = fit(pp_init, h)

# output

TemporalPoissonProcess{Float64}([0.5999999999996618, 1.1400000000005681, 1.7900000000002536])
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