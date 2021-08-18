push!(LOAD_PATH, "../src/")

using Documenter
using PointProcesses

DocMeta.setdocmeta!(
    PointProcesses,
    :DocTestSetup,
    :(using PointProcesses);
    recursive = true,
)

makedocs(;
    modules = [PointProcesses],
    authors = "Guillaume Dalle",
    repo = "https://github.com/gdalle/PointProcesses.jl/blob/{commit}{path}#{line}",
    sitename = "PointProcesses.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://gdalle.github.io/PointProcesses.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => [
            "Event history" => "tuto/history.md",
            "Point processes" => "tuto/point_processes.md",
            "Built-in models" => "tuto/models.md",
            "Markov chains" => "tuto/markov.md",
            "Hidden Markov models" => "tuto/hmm.md",
            "Utilities" => "tuto/utils.md",
        ],
        "API reference" => "api.md",
        "Roadmap" => "roadmap.md",
    ],
)

deploydocs(; repo = "github.com/gdalle/PointProcesses.jl", devbranch = "master")
