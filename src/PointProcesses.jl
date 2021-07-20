"""
A package for point process modeling, simulation and inference.
"""
module PointProcesses

using ComponentArrays
using Distributions
using ForwardDiff
using GalacticOptim
using LinearAlgebra
using Optim
using Plots
using Quadrature
using Random: AbstractRNG, GLOBAL_RNG

using NamedTupleTools: ntfromstruct

# General stuff
include("history.jl")

# Markov processes
include("markov.jl")
include("hmm.jl")
include("forward_backward.jl")

# Point processes
include("point_process.jl")
include("poisson.jl")
include("hawkes.jl")
include("self_correcting.jl")

# Utils
include("utils.jl")

export History
export nb_events, has_events, duration

export PointProcess, MultivariatePointProcess, Parameter
export get_θ
export intensity, mark_distribution, ground_intensity, ground_intensity_bound
export integrated_ground_intensity, logpdf, fit
export rand

export MultivariatePoissonProcess, MultivariateHawkesProcess

export logsumexp, scatter

end
