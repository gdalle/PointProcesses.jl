using Documenter
using PointProcesses
using Test

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
