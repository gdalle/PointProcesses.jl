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
        "Tutorial" => "tutorial.md",
        "Components" => [
            "Event History" => "history.md",
            "Markov models" => [
                "Markov chains" => "markov.md",
                "Hidden Markov models" => "hmm.md",
            ],
            "Point Processes" => [
                "General functions" => "point_processes.md",
                "Built-in models" => "models.md",
            ],
            "Utilities" => "utils.md",
        ],
    ],
)

# deploydocs(; repo = "github.com/gdalle/PointProcesses.jl", devbranch = "main")
