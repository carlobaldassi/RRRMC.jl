module Mixed

using ExtractMacro
using ..Common
using ..Interface

export GraphMixed

import ..Interface: energy, delta_energy, neighbors, update_cache!


struct GraphMixed{ET} <: SimpleGraph{ET}
    N::Int
    graphs::Vector{AbstractGraph}
end

"""
    GraphMixed(graphs::AbstractGraph...)

A `SimpleGraph` mixing 2 or more base graphs into a single model. 
The energy function of a mixed graph is just the sum of the energy functions
of the individual graphs in `graphs`. 

***Usage Example***:

    X = RRRMC.GraphMixed(RRRMC.GraphRRG(10, 3), 
                         RRRMC.GraphSK(10))

    X = RRRMC.GraphMixed(RRRMC.GraphSK(31),
                         RRRMC.GraphPercStep(31, 30),
                         RRRMC.GraphSAT(31, 3, 4.2))
"""
function GraphMixed(graphs::AbstractGraph...)
    length(graphs) >= 2 || throw(ArgumentError("at least 2 base graphs needed"))
    N = getN(graphs[1])
    all(getN(x) == N for x in graphs) || throw(ArgumentError("same N for all graphs required"))
    ET = typeof(sum(zero(energy_type(x)) for x in graphs))
    gs = AbstractGraph[x for x in graphs]
    return GraphMixed{ET}(N, gs)
end

function energy(X::GraphMixed{ET}, C::Config) where ET
    sum(energy(x, C) for x in X.graphs)::ET
end

function delta_energy(X::GraphMixed{ET}, C::Config, move::Int) where ET
    sum(delta_energy(x, C, move) for x in X.graphs)::ET
end

function neighbors(X::GraphMixed, i::Int)
    union(neighbors.(X.graphs, i)...)
end 

function update_cache!(X::GraphMixed, C::Config, move::Int)
    foreach(x->update_cache!(x, C, move), X.graphs)
end 

end #module
