"""
A package for point process modeling, simulation and inference.
"""
module PointProcesses

## Imports

using Distributions
using ForwardDiff
using GalacticOptim
using LinearAlgebra
using LogExpFunctions
using NamedTupleTools
using Optim
using Parameters
using Plots
using Quadrature
using Random
using StatsBase
using StatsPlots
using TransformVariables
using Zygote

using Random: GLOBAL_RNG

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
export DiscreteMarkovChain, transition_matrix
export initial_distribution, stationary_distribution

include("markov/continuous_time.jl")
export ContinuousMarkovChain, rate_matrix, rate_diag, discretize

# Point processes

include("point_processes/abstract.jl")
export AbstractPointProcess, build_transform

include("point_processes/temporal.jl")
export TemporalPointProcess, BoundedTemporalPointProcess

include("point_processes/intensity.jl")
export intensity, mark_distribution, ground_intensity, ground_intensity_bound
export MultivariateTemporalPointProcess, all_marks, all_mark_probabilities

include("point_processes/learning.jl")
export integrated_ground_intensity, check_residuals

include("point_processes/ogata.jl")

# Hidden Markov models

include("hmm/hmm.jl")
export HiddenMarkovModel
export transitions, emissions, emission

include("hmm/forward_backward.jl")
export forward_nolog!, forward_log!
export backward_nolog!, backward_log!
export forward_backward_nolog!, forward_backward_log!

include("hmm/baum_welch.jl")
export update_observation_likelihood!, baum_welch_step!, baum_welch!, baum_welch


# Models

include("models/poisson.jl")
export PoissonProcess

include("models/poisson_inhomogeneous.jl")
export InhomogeneousPoissonProcess

include("models/poisson_multivariate.jl")
export MultivariatePoissonProcess

include("models/poisson_multivariate_naive.jl")
export NaiveMultivariatePoissonProcess

include("models/hawkes.jl")
export MultivariateHawkesProcess

# Utils

include("utils/utils.jl")
export all_minus_inf, all_plus_inf, all_zeros, all_nan

include("utils/categorical.jl")
export fit_map
export CategoricalPrior

include("utils/mvcategorical.jl")
export MvCategorical, MvCategoricalStats, MvCategoricalPrior

include("utils/plot.jl")
export plot_events, plot_intensity, qqplot_interevent_times

include("utils/randvals.jl")
export uniformprobvec, randprobvec
export uniformtransmat, randtransmat

end
