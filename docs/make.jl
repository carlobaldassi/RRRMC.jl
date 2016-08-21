using Documenter, RRRMC

makedocs(
    modules = [RRRMC],
    # format  = Documenter.Formats.HTML,
    # pages   = Any[
    #     "Home" => "index.md",
    #     "algorithms.md",
    #     "graph-types.md",
    #     "graphs-builtin.md",
    #     "interface.md"
    #    ]
    )

deploydocs(
    deps   = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo   = "github.com/carlobaldassi/RRRMC.jl.git",
    julia  = "release"
)

