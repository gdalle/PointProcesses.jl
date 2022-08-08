using Aqua
using Documenter
using Distributions
using JuliaFormatter
using PointProcesses
using Random
using Statistics
using Test

Random.seed!(63)

DocMeta.setdocmeta!(PointProcesses, :DocTestSetup, :(using PointProcesses); recursive=true)

@testset verbose = true "PointProcesses.jl" begin
    @testset verbose = true "Code quality (Aqua.jl)" begin
        Aqua.test_all(PointProcesses; ambiguities=false)
    end
    @testset verbose = true "Formatting" begin
        @test format(PointProcesses; verbose=true, overwrite=false)
    end
    @testset verbose = true "Doctests" begin
        doctest(PointProcesses)
    end
    @testset verbose = true "History" begin
        include("history.jl")
    end
    @testset verbose = true "Poisson" begin
        @testset verbose = true "Multivariate" begin
            include("multivariate_poisson_process.jl")
        end
        @testset verbose = true "Marked" begin
            include("marked_poisson_process.jl")
        end
    end
end
