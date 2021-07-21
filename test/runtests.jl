using Documenter
using PointProcesses
using Test

DocMeta.setdocmeta!(
    PointProcesses,
    :DocTestSetup,
    :(using PointProcesses);
    recursive = true,
)

@testset "PointProcesses.jl doctests" begin
    doctest(PointProcesses)
end
