using Documenter
using Distributions
using PointProcesses
using Quadrature
using Random
using Statistics
using Test

Random.seed!(96)

@testset verbose = true "PointProcesses.jl" begin
    include("temporal.jl")
    include("models.jl")
    include("markov.jl")
    include("hmm.jl")
    # include("mmpp.jl")
end
