using Documenter
using PointProcesses
using Random
using Test

Random.seed!(63)

DocMeta.setdocmeta!(
    PointProcesses,
    :DocTestSetup,
    :(using PointProcesses);
    recursive = true,
)

@testset "Doctests" begin
    doctest(PointProcesses)
end

include("poisson.jl")
