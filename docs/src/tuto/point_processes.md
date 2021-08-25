# General framework for point processes

```@setup point_processes
using ForwardDiff, GalacticOptim, Optim, Random
Random.seed!(63)
```

The main goal of this package lies in point process simulation and inference. All point processes are subtypes of [`AbstractPointProcess{L,M}`](@ref), where `L` is the type of event locations and `M` is the type of event marks.

## Temporal point processes

Temporal point processes have specific features which make it possible to define generic functions to simulate or estimate them.[^Rasmussen_2018]

[^Rasmussen_2018]: Rasmussen, J. G. (2018), “Lecture Notes: Temporal Point Processes and the Conditional Intensity Function,” arXiv:1806.00221 [stat].


To implement your own process, you only have to define a subtype of [`TemporalPointProcess`](@ref) and implement the necessary methods: [`intensity`](@ref), [`mark_distribution`](@ref), [`ground_intensity`](@ref) and [`ground_intensity_bound`](@ref).
As long as these methods exist, the default simulation and inference routines should work, but they can be made much more efficient using custom implementations (for instance to avoid numerical optimization or gradient optimization).

As an example, we included a naive implementation of a Poisson process with categorical mark distribution, called [`NaiveMultivariatePoissonProcess`](@ref). Looking at its source code may help you understand the requirements of the interface. In the meantime, we can check that the implementation works:

```@repl point_processes
using GalacticOptim, PointProcesses, Random; Random.seed!(63)

pp = NaiveMultivariatePoissonProcess([1., 2., 3.])
h = rand(pp, 0., 100.)
pp_init = NaiveMultivariatePoissonProcess(ones(3))
pp_est1 = fit(pp_init, h)
intensity(pp_est1)
```

We can also tune the automatic differentiation method and the optimization algorithm:

```@repl point_processes
using ForwardDiff, Optim

pp_est2 = fit(pp_init, h, adtype=GalacticOptim.AutoForwardDiff(), alg=Optim.LBFGS())
intensity(pp_est2)
```