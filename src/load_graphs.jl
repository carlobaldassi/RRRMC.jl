# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

const graphs_dir = joinpath(@__DIR__, "graphs")
## load all julia files in the graphs_dir, without descending in subdirectories
## silently skip files starting with an underscore
const valid_graph_name = r"^([^/_][^/]*)\.jl$"

macro include_graph(filename)
    if !occursin(valid_graph_name, filename)
        startswith(filename, '_') || @warn("Unrecogniezd file $filename, skipping")
        return :()
    end
    modname = Symbol(replace(filename, valid_graph_name => s"\1"))
    quote
        include(joinpath(graphs_dir, $(esc(filename))))
        using .$modname
    end
end

@include_graph "QT.jl"
@include_graph "RE.jl"
@include_graph "LE.jl"
@include_graph "TLE.jl"
@include_graph "AddFields.jl"

for filename in filter(f->endswith(f, ".jl") && f ∉ ["QT.jl", "RE.jl", "LE.jl", "TLE.jl", "AddFields.jl"], readdir(graphs_dir))
    @eval @include_graph $filename
end
