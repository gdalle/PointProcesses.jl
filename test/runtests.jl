using Documenter
using ForwardDiff
using GalacticOptim
using MeasureTheory
using Optim
using PointProcesses
using Random
using Test

Random.seed!(63)

include("history.jl")

include("markov.jl")

include("point_processes.jl")

include("models.jl")

include("hmm.jl")

include("utils.jl")

include("doctests.jl")
