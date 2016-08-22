using Documenter, RRRMC

makedocs(
    modules  = [RRRMC],
    # format   = Documenter.Formats.HTML,
    # sitename = "RRRMC.jl",
    # pages    = Any[
    #     "Home" => "index.md",
    #     "Sampling algorithms" => "algorithms.md",
    #     "Graph types" => "graph-types.md",
    #     "Built-in graphs" => "graphs-builtin.md",
    #     "Graphs interface" => "interface.md"
    #    ]
    )

deploydocs(
    deps   = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo   = "github.com/carlobaldassi/RRRMC.jl.git",
    julia  = "release"
)

