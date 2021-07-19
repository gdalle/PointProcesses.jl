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

using NamedTupleTools: ntfromstruct

include("history.jl")
include("point_process.jl")
include("utils.jl")
include("ogata.jl")
include("poisson.jl")

export History
export nb_events, has_events, duration

export PointProcess
export get_θ
export intensity, mark_distribution, ground_intensity, ground_intensity_bound

export PoissonProcess

export integrated_ground_intensity
export rand, logpdf, fit

export scatter

end
