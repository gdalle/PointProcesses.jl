push!(LOAD_PATH, "../src/")

using Documenter
using DocumenterCitations
using PointProcesses

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:authoryear)

DocMeta.setdocmeta!(PointProcesses, :DocTestSetup, :(using PointProcesses); recursive=true)

makedocs(
    bib;
    modules=[PointProcesses],
    authors="Guillaume Dalle",
    repo="https://github.com/gdalle/PointProcesses.jl/blob/{commit}{path}#{line}",
    sitename="PointProcesses.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://gdalle.github.io/PointProcesses.jl",
        assets=String[],
    ),
    pages=["Home" => "index.md", "API reference" => "api.md"],
)

deploydocs(; repo="github.com/gdalle/PointProcesses.jl", devbranch="main")
