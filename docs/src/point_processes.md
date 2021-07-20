# Point processes

## General structure

```@docs
PointProcess
Parameter
get_θ
```

## Implementing new processes

```@docs
intensity
mark_distribution
ground_intensity
ground_intensity_bound
```

## Built-in models

```@docs
MultivariatePoissonProcess
MultivariateHawkesProcess
```

## Simulation

```@docs
rand(::PointProcess, ::Real, ::Real)
```

## Learning

```@docs
integrated_ground_intensity
logpdf
fit
```