using Documenter
using PointProcesses
using Random
using Test

@testset verbose = true "PointProcesses.jl" begin
    include("history.jl")
    include("markov.jl")
    include("point_processes.jl")
    include("models.jl")
    # include("hmm.jl")
    # include("utils.jl")
    # include("doctests.jl")
end
