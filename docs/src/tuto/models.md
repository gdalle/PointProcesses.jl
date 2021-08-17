# Built-in models

```@meta
DocTestSetup = quote
    using PointProcesses
end
```

```jldoctest models
julia> using Random

julia> Random.seed!(63);
```

## Poisson processes

We provide a number of built-in models, starting with basic Poisson processes on the real line.

### Multivariate Poisson processes

```jldoctest models
julia> pp = MultivariatePoissonProcess(λ = [0.5, 1., 2.]);

julia> history = rand(pp, 0., 1000.);

julia> pp_est = fit(MultivariatePoissonProcess, history);

julia> round.(pp_est.λ, digits=2)
3-element Vector{Float64}:
 0.49
 0.96
 2.0
```

### General Poisson processes

> Todo