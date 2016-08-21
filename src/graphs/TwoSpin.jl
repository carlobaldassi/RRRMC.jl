module TwoSpin

using ExtractMacro
using ..Interface

export GraphTwoSpin

import ..Interface: energy, delta_energy, neighbors, allΔE

type GraphTwoSpin <: DiscrGraph{Int}
end

@doc """
    GraphTwoSpin() <: DiscrGraph

A trivial `DiscrGraph` type with 2 spins inteacting ferromagnetically
(\$J=1\$), without fields.

Only useful for testing/debugging purposes.
""" -> GraphTwoSpin()

getN(::GraphTwoSpin) = 2

function energy(X::GraphTwoSpin, C::Config)
    @extract C : N s
    @assert N == 2
    return -(2s[1] - 1) * (2s[2] - 1)
end

function delta_energy(X::GraphTwoSpin, C::Config, move::Int)
    @extract C : N s
    @assert 1 ≤ move ≤ N
    @assert N == 2

    return 2 * (2s[1] - 1) * (2s[2] - 1)
end

neighbors(X::GraphTwoSpin, i::Int) = return (3 - i,)
allΔE(::Type{GraphTwoSpin}) = (2,)

end
