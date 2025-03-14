using Aqua
using Documenter
using Distributions
using JuliaFormatter
using PointProcesses
using Random
using Statistics
using LinearAlgebra
using Test

Random.seed!(63)

DocMeta.setdocmeta!(PointProcesses, :DocTestSetup, :(using PointProcesses); recursive=true)

@testset verbose = true "PointProcesses.jl" begin
    @testset verbose = false "Code quality (Aqua.jl)" begin
        Aqua.test_all(PointProcesses; ambiguities=false, deps_compat=(; check_extras=false))
    end
    @testset verbose = false "Formatting" begin
        @test format(PointProcesses; verbose=false, overwrite=false)
    end
    @testset verbose = false "Doctests" begin
        doctest(PointProcesses)
    end
    @testset verbose = true "History" begin
        include("history.jl")
    end
    @testset verbose = true "Bounded" begin
        include("bounded_point_process.jl")
    end
    @testset verbose = true "Poisson" begin
        @testset verbose = true "Multivariate" begin
            include("multivariate_poisson_process.jl")
        end
        @testset verbose = true "Marked" begin
            include("marked_poisson_process.jl")
        end
        @testset verbose = true "Multivariate as marked" begin
            include("multivariate_as_marked.jl")
        end
    end
end
