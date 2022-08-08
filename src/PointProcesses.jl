"""
A package for temporal point process modeling, simulation and inference.
"""
module PointProcesses

# Imports

using DensityInterface: DensityInterface, densityof, logdensityof
using Distributions: Distributions, UnivariateDistribution, MultivariateDistribution
using Distributions: Categorical, Exponential, Poisson, Uniform
using Distributions: fit_mle, suffstats
using LinearAlgebra
using Random
using Random: GLOBAL_RNG

## Hidden names

# Exports

## Reexports

export rand  # Base
export fit_mle, suffstats  # Distributions
export logdensityof, densityof # DensityInterface
export fit_map

## History

export History
export event_times, event_marks, min_time, max_time
export nb_events, has_events, duration
export time_change, split_into_chunks

## Point processes

export AbstractPointProcess
export BoundedPointProcess
export ground_intensity, mark_distribution
export intensity, log_intensity
export ground_intensity_bound
export integrated_ground_intensity
export check_residuals

## Models

export AbstractPoissonProcess
export MultivariatePoissonProcess
export MarkedPoissonProcess

# Includes

include("history.jl")
include("abstract_point_process.jl")
include("simulation.jl")
include("bounded.jl")

include("poisson/abstract_poisson_process.jl")
include("poisson/simulation.jl")

include("poisson/multivariate/multivariate_poisson_process.jl")
include("poisson/multivariate/suffstats.jl")
include("poisson/multivariate/prior.jl")
include("poisson/multivariate/fit.jl")

include("poisson/marked/marked_poisson_process.jl")
include("poisson/marked/suffstats.jl")
include("poisson/marked/fit.jl")

end
