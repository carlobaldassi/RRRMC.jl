using Documenter, RRRMC

CIbuild = get(ENV, "CI", nothing) == "true"

makedocs(
    modules  = [RRRMC],
    format   = Documenter.HTML(prettyurls = CIbuild),
    sitename = "RRRMC.jl",
    pages    = Any[
        "Home" => "index.md",
        "Sampling algorithms" => "algorithms.md",
        "Graph types" => "graph-types.md",
        "Built-in graphs" => "graphs-builtin.md",
        "Graphs interface" => "interface.md"
       ]
    )

deploydocs(
    repo   = "github.com/carlobaldassi/RRRMC.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
