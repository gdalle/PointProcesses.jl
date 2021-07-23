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

## Includes

# General stuff
include("history.jl" )

# Markov processes
include("markov/matrices.jl")
include("markov/markov.jl")
include("markov/markov_discrete.jl")
include("markov/markov_continuous.jl")

# hidden Markov models
include("hmm/hmm.jl")
include("hmm/forward_backward.jl")
include("hmm/baum_welch.jl")

# Point processes
include("pp/point_process.jl")

# Temporal point processes
include("tpp/temporal_point_process.jl")
include("tpp/multivariate.jl")
include("tpp/intensity.jl")
include("tpp/ogata.jl")
include("tpp/learning.jl")
include("tpp/poisson.jl")
include("tpp/hawkes.jl")

# Utils
include("utils/utils.jl")
include("utils/plot.jl")

## Exports

export rand, logpdf, fit, scatter

# .

export History
export nb_events, has_events, duration, time_change

# markov

export AbstractMarkovChain, DiscreteMarkovChain, ContinuousMarkovChain
export nstates, stationary_distribution

# hmm

export HiddenMarkovModel
export forward_log!, forward_nolog!, backward_log!, backward_nolog!

# pp

export PointProcess
export params, build_transform

# tpp

export TemporalPointProcess, Bounded
export Parameter, params
export MultivariatePointProcess, all_marks
export intensity, mark_distribution, ground_intensity, ground_intensity_bound
export integrated_ground_intensity
export check_residuals

export PoissonProcess, HawkesProcess

# utils

export plot_events, plot_intensity, qqplot_interevent_times
export logsumexp

end
