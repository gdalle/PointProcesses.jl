using Documenter
using DocumenterCitations
using PointProcesses

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:authoryear)

DocMeta.setdocmeta!(PointProcesses, :DocTestSetup, :(using PointProcesses); recursive=true)

makedocs(;
    modules=[PointProcesses],
    authors="Guillaume Dalle",
    sitename="PointProcesses.jl",
    format=Documenter.HTML(),
    pages=["Home" => "index.md", "API reference" => "api.md"],
    plugins=[bib],
)

deploydocs(; repo="github.com/gdalle/PointProcesses.jl", devbranch="main")
