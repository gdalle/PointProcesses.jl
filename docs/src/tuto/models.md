# Built-in models

```@meta
DocTestSetup = quote
    using PointProcesses
end
```

## Poisson processes

### Multivariate Poisson processes

Multivariate Poisson processes are a superposition of a finite number of univariate Poisson processes with consecutive integer marks.

```jldoctest
julia> using Random; Random.seed!(63);

julia> pp = MultivariatePoissonProcess([0.5, 1., 2.]);

julia> history = rand(pp, 0., 1000.);

julia> pp_est = fit(MultivariatePoissonProcess, history);

julia> round.(pp_est.λ, digits=2)
3-element Vector{Float64}:
 0.49
 0.96
 2.0
```

### General Poisson processes

We also provide a more general [`PoissonProcess`](@ref) structure with a single intensity value and an arbitrary mark distribution. For instance, the `Categorical` distribution enables us to recreate the multivariate Poisson process from the previous section:

```jldoctest
julia> using Random, Distributions; Random.seed!(63);

julia> mark_dist = Categorical([0.5, 1., 2.] / 3.5);

julia> pp = PoissonProcess{Int}(3.5, mark_dist);

julia> h = rand(pp, 0., 100.);

julia> pp_est = fit(PoissonProcess{Int,Categorical}, h);

julia> round.(probs(mark_distribution(pp_est)), digits=2)
3-element Vector{Float64}:
 0.14
 0.28
 0.57
```