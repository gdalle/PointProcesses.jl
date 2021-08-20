# General framework for point processes

```@meta
DocTestSetup = quote
    using PointProcesses
end
```

The main goal of this package lies in point process simulation and inference. All point processes are subtypes of [`AbstractPointProcess{L,M}`](@ref), where `L` is the type of event locations and `M` is the type of event marks.

## Temporal point processes

Temporal point processes have specific features which make it possible to define generic functions to simulate or estimate them.

To implement your own process, you only have to define a subtype of [`TemporalPointProcess`](@ref) and implement the necessary methods: [`intensity`](@ref), [`mark_distribution`](@ref), [`ground_intensity`](@ref) and [`ground_intensity_bound`](@ref).
As long as these methods exist, the default simulation and inference routines should work, but they can be made much more efficient using custom implementations (for instance to avoid numerical optimization or gradient optimization).

As an example, we included a naive implementation of a Poisson process with categorical mark distribution, called [`NaiveMultivariatePoissonProcess`](@ref). Looking at its source may help you understand the requirements of the interface.