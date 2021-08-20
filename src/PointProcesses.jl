"""
A package for point process modeling, simulation and inference.
"""
module PointProcesses

## Imports

using DataStructures
using ForwardDiff
using GalacticOptim
using LinearAlgebra
using LogExpFunctions
using MeasureTheory
using NamedTupleTools
using OffsetArrays
using Optim
using Plots
using Quadrature
using Random
using SimpleTraits
using StatsPlots
using TransformVariables
using Zygote

# Hidden names
using Random: GLOBAL_RNG

# Functions to extend
import Base: eltype, length, rand
import Distributions: fit, fit_mle, logpdf, pdf, suffstats
import MeasureTheory: density, logdensity, sampletype, testvalue

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
export TemporalHistory
export event_times, event_marks, min_time, max_time
export nb_events, has_events, duration, time_change

# Markov processes

include("markov/abstract.jl")
export AbstractMarkovChain
export nb_states
export AbstractMarkovChainPrior

include("markov/discrete_time.jl")
export DiscreteMarkovChain
export initial_distribution, transition_matrix, stationary_distribution
export DiscreteMarkovChainPrior, DiscreteMarkovChainStats

include("markov/continuous_time.jl")
export ContinuousMarkovChain
export rate_matrix, rate_diag, discretize
export ContinuousMarkovChainPrior, ContinuousMarkovChainStats

# Point processes

include("point_processes/abstract.jl")
export AbstractPointProcess, build_transform

include("point_processes/temporal.jl")
export TemporalPointProcess, BoundedTemporalPointProcess
export intensity, mark_distribution, ground_intensity, ground_intensity_bound
export MultivariateTemporalPointProcess, all_marks, all_mark_probabilities
export integrated_ground_intensity, check_residuals

# Models

include("models/poisson.jl")
export PoissonProcess

include("models/poisson_multivariate_naive.jl")
export NaiveMultivariatePoissonProcess

# Hidden Markov models

include("hmm/hmm.jl")
export HiddenMarkovModel
export transitions, emissions, emission

include("hmm/baum_welch.jl")
export forward_nolog!, forward_log!
export backward_nolog!, backward_log!
export forward_backward_nolog!, forward_backward_log!
export update_obs_pdf!, update_obs_logpdf!
export baum_welch!, baum_welch

include("hmm/mmpp.jl")
export MarkovModulatedPoissonProcess

include("hmm/ryden.jl")
export forward_backward, ryden

# Utils

include("utils/overflow.jl")
export all_minus_inf, all_plus_inf, all_zeros, all_nan

include("utils/dists.jl")

include("utils/categorical.jl")

include("utils/plot.jl")
export plot_events, plot_intensity, qqplot_interevent_times

include("utils/randvals.jl")
export uniformprobvec, randprobvec
export uniformtransmat, randtransmat

end
