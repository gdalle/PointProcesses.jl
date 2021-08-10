"""
A package for point process modeling, simulation and inference.
"""
module PointProcesses

## Imports

using Distributions
using ForwardDiff
using LinearAlgebra
using Optim
using Plots
using Quadrature
using Random
using StatsBase
using StatsPlots
using TransformVariables

using Random: AbstractRNG, GLOBAL_RNG
using NamedTupleTools: ntfromstruct

## Includes and exports

# Reexports

export eltype, rand  # Base
export params  # StatsBase
export logpdf, fit, fit_mle, suffstats  # Distributions
export plot, scatter  # Plots

# History

include("history/abstract.jl")
export AbstractHistory

include("history/temporal.jl")
export TemporalHistory, nb_events, has_events, duration, time_change

# Markov processes

include("markov/abstract.jl")
export AbstractMarkovChain, nstates

include("markov/discrete_time.jl")
export DiscreteMarkovChain, stationary_distribution

include("markov/continuous_time.jl")
export ContinuousMarkovChain, discretize

# Hidden Markov models

include("hmm/hmm.jl")
export HiddenMarkovModel, update_observation_likelihood!

include("hmm/forward_backward.jl")
export forward_nolog!, forward_log!
export backward_nolog!, backward_log!
export forward_backward_nolog!, forward_backward_log!

include("hmm/baum_welch.jl")
export expectation_maximization_step!

# Point processes

include("point_processes/abstract.jl")
export AbstractPointProcess

include("point_processes/temporal.jl")
export TemporalPointProcess, BoundedTemporalPointProcess

include("point_processes/intensity.jl")
export intensity, mark_distribution, ground_intensity, ground_intensity_bound

include("point_processes/temporal_multivariate.jl")
export MultivariateTemporalPointProcess, all_marks, all_mark_probabilities

include("point_processes/learning.jl")
export integrated_ground_intensity, check_residuals

include("point_processes/ogata.jl")

# Models

include("models/temporal_poisson.jl")
export TemporalPoissonProcess, build_transform

include("models/temporal_hawkes.jl")
export TemporalHawkesProcess

# Utils

include("utils/utils.jl")
export logsumexp

include("utils/plot.jl")
export plot_events, plot_intensity, qqplot_interevent_times

end
