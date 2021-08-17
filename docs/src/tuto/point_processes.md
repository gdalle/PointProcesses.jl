# General framework for point processes

```@meta
DocTestSetup = quote
    using PointProcesses
end
```

## Working with temporal point processes

The main goal of this package lies in point process simulation and inference. All point processes are subtypes of [`AbstractPointProcess{L,M}`](@ref), where `L` is the type of event locations and `M` is the type of event marks.

To implement your own process, you only have to define a subtype of [`TemporalPointProcess`](@ref) and write the necessary methods: [`intensity`](@ref), [`mark_distribution`](@ref), [`ground_intensity`](@ref) and [`ground_intensity_bound`](@ref).

As long as these methods exist, the default simulation and inference routines should work, but they can be made much more efficient using custom implementations.

As an example, we included a naive re-implementation of 
[`MultivariatePoissonProcess`](@ref), called [`NaiveMultivariatePoissonProcess`](@ref). Looking at its source may help you understand the requirements of the interface.