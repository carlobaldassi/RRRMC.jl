# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

const graphs_dir = joinpath(dirname(@__FILE__), "graphs")
const valid_graph_name = r"^([^/]+)\.jl$"

macro include_graph(filename)
    ismatch(valid_graph_name, filename) || (warn("Unrecogniezd file $filename, skipping"); return :())
    modname = Symbol(replace(filename, valid_graph_name, s"\1"))
    quote
        include(joinpath(graphs_dir, $(esc(filename))))
        using .$(esc(modname))
    end
end

@include_graph "QT.jl"
@include_graph "RE.jl"
@include_graph "LE.jl"

for filename in filter(f->endswith(f, ".jl") && f ∉ ["QT.jl", "RE.jl", "LE.jl"], readdir(graphs_dir))
    @eval @include_graph $filename
end
