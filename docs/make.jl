using Documenter, RRRMC

makedocs(
    modules  = [RRRMC],
    format   = :html,
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
    julia  = "0.5"
)
