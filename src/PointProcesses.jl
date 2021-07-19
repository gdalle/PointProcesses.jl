module PointProcesses

using ComponentArrays
using Distributions
using ForwardDiff
using GalacticOptim
using Optim
using ParameterHandling
using Plots
using Quadrature
using StatsFuns

include("history.jl")
include("point_process.jl")
include("utils.jl")
include("ogata.jl")
include("poisson.jl")
include("learning.jl")

export History
export nb_events, has_events

export PointProcess
export default_params
export intensity, mark_distribution, ground_intensity, ground_intensity_bound

export PoissonProcess

export integrated_ground_intensity
export rand, logpdf, fit

export scatter

end
