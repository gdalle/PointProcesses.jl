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
           emissions = [Normal(1, 0.3), Normal(2, 0.3)]
       )
HiddenMarkovModel{DiscreteMarkovChain{Float64}, Normal{Float64}}
  transitions: DiscreteMarkovChain{Float64}
  emissions: Array{Normal{Float64}}((2,))

julia> states, observations = rand(hmm, 1000);

julia> hmm_init = HiddenMarkovModel(
           transitions = DiscreteMarkovChain(π0 = randprobvec(2), P = randtransmat(2)),
           emissions = [Normal(rand(), 1), Normal(rand(), 1)]
       );

julia> hmm_est, logL_evolution = baum_welch(hmm_init, observations, iter=100);

julia> minimum(diff(logL_evolution))
-8.185452315956354e-12

julia> round.(transition_matrix(hmm_est), digits=2)
2×2 Matrix{Float64}:
 0.83  0.17
 0.11  0.89

julia> emissions(hmm_est)
2-element Vector{Normal{Float64}}:
 Distributions.Normal{Float64}(μ=2.007774167238364, σ=0.28677264095786803)
 Distributions.Normal{Float64}(μ=0.9764685756929417, σ=0.30552164217241995)
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
 2.0
```

### General Poisson processes

> Todo

### Implementing your own models

To implement your own process, you only have to define a subtype of [`TemporalPointProcess`](@ref) and write the necessary methods: [`intensity`](@ref), [`mark_distribution`](@ref), [`ground_intensity`](@ref) and [`ground_intensity_bound`](@ref).

As long as these methods exist, the default simulation and inference routines should work, but they can be made much more efficient using custom implementations.

As an example, we included a naive re-implementation of 
[`MultivariatePoissonProcess`](@ref), called [`NaiveMultivariatePoissonProcess`](@ref). Looking at its source may help you understand the requirements of the interface.