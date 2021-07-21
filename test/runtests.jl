using Documenter
using PointProcesses
using Test

DocMeta.setdocmeta!(
    PointProcesses,
    :DocTestSetup,
    :(using PointProcesses; using Random);
    recursive = true,
)

@testset "PointProcesses.jl doctests" begin
    doctest(PointProcesses)
end
