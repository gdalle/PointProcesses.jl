"""
A package for point process modeling, simulation and inference.
"""
module PointProcesses

# Imports

using DataStructures
using DensityInterface
using Distributions
using FillArrays
using LinearAlgebra
using LogExpFunctions
using OffsetArrays
using ProgressMeter
using Random
using Random: GLOBAL_RNG
using UnicodePlots
using UnPack

## Hidden names

# Exports

## Reexports

export rand  # Base
export fit, fit_mle, suffstats  # Distributions
export logdensityof # DensityInterface

export fit_map

## History

export History
export event_times, event_marks, min_time, max_time, times_and_marks
export nb_events, has_events, duration
export time_change, split_into_chunks

## Markov processes

export DiscreteMarkovChain, DiscreteMarkovChainPrior
export initial_distribution, transition_matrix, stationary_distribution
export nb_states

export ContinuousMarkovChain, ContinuousMarkovChainPrior
export rate_matrix, rate_diag, discretize_chain

## Point processes

export TemporalPointProcess, BoundedTemporalPointProcess
export intensity, log_intensity, mark_distribution
export ground_intensity, ground_intensity_bound
export integrated_ground_intensity
export check_residuals

## Models

export PoissonProcess
export MultivariatePoissonProcess

## Hidden Markov models

export HiddenMarkovModel
export transitions, emissions, emission

export forward!, backward!, update_obs_density!
export forward_log!, backward_log!, update_obs_logdensity!
export baum_welch_multiple_sequences, baum_welch_multiple_sequences_log, baum_welch

# export MarkovModulatedPoissonProcess

# export forward_backward, ryden

## Utils

export uniformprobvec, randprobvec
export uniformtransmat, randtransmat

# Includes

## Unconditional

include("temporal/history.jl")
include("temporal/point_processes.jl")
include("temporal/simulation.jl")

include("markov/discrete_time.jl")
include("markov/continuous_time.jl")

include("models/poisson_generic.jl")
include("models/poisson_multivariate.jl")

include("hmm/hmm.jl")
include("hmm/forward_backward.jl")
include("hmm/forward_backward_log.jl")
include("hmm/baum_welch.jl")

# include("mmpp/mmpp.jl")
# include("mmpp/ryden.jl")

include("utils/overflow.jl")
include("utils/randvals.jl")
include("utils/plot.jl")
include("utils/product_dist.jl")

end
