"""
A package for point process modeling, simulation and inference.
"""
module PointProcesses

using Base: throw_setindex_mismatch
using ComponentArrays
using Distributions
using ForwardDiff
using GalacticOptim
using LinearAlgebra
using Optim
using Plots
using Quadrature
using Random
using StatsBase

using Random: AbstractRNG, GLOBAL_RNG
using NamedTupleTools: ntfromstruct

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
include("pp/multivariate.jl")
include("pp/poisson.jl")
include("pp/hawkes.jl")
include("pp/intensity.jl")
include("pp/ogata.jl")
include("pp/learning.jl")

# Utils
include("utils/utils.jl")

export rand, logpdf, fit, scatter

export History
export nb_events, has_events, duration

export AbstractMarkovChain, DiscreteMarkovChain, ContinuousMarkovChain
export nstates, stationary_distribution

export HiddenMarkovModel
export forward_log!, forward_nolog!, backward_log!, backward_nolog!

export PointProcess, TimedPointProcess
export Parameter, params
export MultivariatePointProcess, all_marks
export intensity, mark_distribution, ground_intensity, ground_intensity_bound
export integrated_ground_intensity

export MultivariatePoissonProcess, MultivariateHawkesProcess

export logsumexp

end
