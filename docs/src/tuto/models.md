# Built-in models

## Poisson processes

We provide a general [`PoissonProcess`](@ref) structure with a single intensity value and an arbitrary mark distribution.

```@repl models
using MeasureTheory, PointProcesses, Random; Random.seed!(63)

mark_dist = Dists.Categorical([0.1, 0.3, 0.6])
pp = PoissonProcess(5., mark_dist)
h = rand(pp, 0., 100.)
pp_est = fit(PoissonProcess{Dists.Categorical}, h)
intensity(pp_est)
Dists.probs(mark_distribution(pp_est))
```