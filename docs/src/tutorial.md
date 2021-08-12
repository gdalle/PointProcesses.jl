# Tutorial

```@meta
DocTestSetup = quote
    using PointProcesses
end
```

In this tutorial we demonstrate the main features of PointProcesses.jl.

```jldoctest tuto
julia> using Random

julia> using Distributions

julia> Random.seed!(63);
```

## Working with event histories

To analyze point processes, we need a way to store their realizations. This is the purpose of the [`AbstractHistory`](@ref) subtypes, of which only [`TemporalHistory`](@ref) is implemented for now.

```jldoctest tuto
julia> history = TemporalHistory([0.2, 0.8, 1.1], ["a", "b", "c"], 0.0, 2.0);

julia> duration(history)
2.0

julia> nb_events(history)
3

julia> nb_events(history, 1.0, 2.0)
1

julia> has_events(history)
true

julia> has_events(history, 1.5, 2.0)
false

julia> push!(history, 1.7, "d")

julia> has_events(history, 1.5, 2.0)
true
```

## Working with Markov processes

Some point processes are based on underlying Markov processes, which is why we provide a basic implementation for them.

### Discrete time Markov chains

```jldoctest tuto
julia> dmc = DiscreteMarkovChain(π0 = [0.3, 0.7], P = [0.9 0.1; 0.2 0.8]);

julia> states = rand(dmc, 1000);

julia> dmc_est = fit(DiscreteMarkovChain, states);

julia> round.(dmc_est.P, digits=1)
2×2 Matrix{Float64}:
 0.9  0.1
 0.2  0.8
```

### Continuous time Markov chains

```jldoctest tuto
julia> cmc = ContinuousMarkovChain(π0 = [0.3, 0.7], Q = [-1. 1.; 2. -2.]);

julia> history = rand(cmc, 0., 1000.);

julia> cmc_est = fit(ContinuousMarkovChain, history);

julia> round.(cmc_est.Q, digits=1)
2×2 Matrix{Float64}:
 -1.0   1.0
  1.9  -1.9
```

### Hidden Markov models

```jldoctest tuto
julia> hmm = HiddenMarkovModel(
           transitions = dmc,
           emissions = [Normal(1, 1), Normal(2, 1)]
       )
HiddenMarkovModel{DiscreteMarkovChain{Float64}, Normal{Float64}}
  transitions: DiscreteMarkovChain{Float64}
  emissions: Array{Normal{Float64}}((2,))

julia> states, observations = rand(hmm, 1000);

julia> hmm_init = HiddenMarkovModel(
           transitions = DiscreteMarkovChain(π0 = randprobvec(2), P = randtransmat(2)),
           emissions = [Normal(rand(), rand()), Normal(rand(), rand())]
       );

julia> hmm_est, logL_evolution = baum_welch(hmm_init, observations, iter=100);

julia> minimum(diff(logL_evolution)) > 0
true

julia> transition_matrix(hmm_est)
2×2 Matrix{Float64}:
 0.199261   0.800739
 0.0336893  0.966311

julia> emissions(hmm_est)
2-element Vector{Normal{Float64}}:
 Normal{Float64}(μ=1.562958876776265, σ=0.2512621023615442)
 Normal{Float64}(μ=1.3325788877249478, σ=1.1437919911321)
```

## Working with temporal point processes

We finally demonstrate the main goal of the package: point process simulation and inference. All point processes are subtypes of [`AbstractPointProcess{L,M}`](@ref), where `L` is the type of event locations and `M` is the type of event marks.

We provide a number of built-in models, starting with basic Poisson processes on the real line.

### Multivariate Poisson processes

```jldoctest tuto
julia> pp = MultivariatePoissonProcess(λ = [0.5, 1., 2.]);

julia> history = rand(pp, 0., 1000.);

julia> pp_est = fit(MultivariatePoissonProcess, history);

julia> round.(pp_est.λ, digits=1)
3-element Vector{Float64}:
 0.5
 1.0
 2.1
```

### General Poisson processes

> Todo

### Implementing your own models

To implement your own process, you only have to define a subtype of [`TemporalPointProcess`](@ref) and write the necessary methods: [`intensity`](@ref), [`mark_distribution`](@ref), [`ground_intensity`](@ref) and [`ground_intensity_bound`](@ref).

As long as these methods exist, the default simulation and inference routines should work, but they can be made much more efficient using custom implementations.

As an example, we included a naive re-implementation of 
[`MultivariatePoissonProcess`](@ref), called [`NaiveMultivariatePoissonProcess`](@ref). Looking at its source may help you understand the requirements of the interface.