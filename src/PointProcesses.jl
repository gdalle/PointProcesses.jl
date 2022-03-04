"""
A package for point process modeling, simulation and inference.
"""
module PointProcesses

# Imports

using DataStructures
using DensityInterface
using Distributions
using LinearAlgebra
using LogExpFunctions
using OffsetArrays
using Random
using Random: GLOBAL_RNG
using Requires
using UnPack

## Hidden names

# Exports

## Reexports

export eltype, rand  # Base
export mean  # Statistics
export fit, fit_mle, suffstats  # Distributions
export plot, scatter  # Plot

## History

export History
export event_times, event_marks, min_time, max_time
export nb_events, has_events, duration
export time_change

## Markov processes

export DiscreteMarkovChain, DiscreteMarkovChainPrior
export initial_distribution, transition_matrix, stationary_distribution

export ContinuousMarkovChain, ContinuousMarkovChainPrior
export rate_matrix, rate_diag, discretize_chain

## Point processes

export TemporalPointProcess, BoundedTemporalPointProcess
export intensity, log_intensity, mark_distribution
export ground_intensity, ground_intensity_bound

## Models

export PoissonProcess
export MultivariatePoissonProcess

## Hidden Markov models

# export HiddenMarkovModel
# export transitions, emissions, emission

# export forward_nolog!, forward_log!
# export backward_nolog!, backward_log!
# export forward_backward_nolog!, forward_backward_log!
# export update_obs_pdf!, update_obs_logpdf!
# export baum_welch!, baum_welch

# export MarkovModulatedPoissonProcess

# export forward_backward, ryden

## Utils

# export all_minus_inf, all_plus_inf, all_zeros, all_nan
# export uniformprobvec, randprobvec
# export uniformtransmat, randtransmat

# Includes

## Unconditional

include("temporal/history.jl")
include("temporal/point_processes.jl")
include("temporal/simulation.jl")
include("temporal/integration.jl")

include("markov/discrete_time.jl")
include("markov/continuous_time.jl")

include("models/poisson_generic.jl")
include("models/poisson_multivariate.jl")

# include("hmm/hmm.jl")
# include("hmm/baum_welch.jl")
# include("hmm/mmpp.jl")
# include("hmm/ryden.jl")

# include("utils/overflow.jl")
# include("utils/categorical.jl")
# include("utils/randvals.jl")

## Conditional dependencies

function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("temporal/plot.jl")
    @require Quadrature = "67601950-bd08-11e9-3c89-fd23fb4432d2" include("temporal/integration.jl")
end

end
