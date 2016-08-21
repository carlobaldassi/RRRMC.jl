module ThreeSpin

using ExtractMacro
using ..Interface

export GraphThreeSpin

import ..Interface: energy, delta_energy, neighbors, allΔE, getN

type GraphThreeSpin <: DiscrGraph{Int}
end

@doc """
    GraphThreeSpin() <: DiscrGraph

A trivial `DiscrGraph` type with 3 spins, ferromagnetic interactions
(\$J=1\$), no fields, and periodic boundary conditions.

Only useful for testing/debugging purposes.
""" -> GraphThreeSpin()

getN(::GraphThreeSpin) = 3

function energy(X::GraphThreeSpin, C::Config)
    @extract C : N s
    @assert N == 3
    return -(2s[1] - 1) * (2s[2] - 1) - (2s[2] - 1) * (2s[3] - 1) - (2s[3] - 1) * (2s[1] - 1)
end

function delta_energy(X::GraphThreeSpin, C::Config, move::Int)
    @extract C : N s
    @assert 1 ≤ move ≤ N
    @assert N == 3

    ΔE = 0
    1 ≤ move ≤ 2 && (ΔE += 2 * (2s[1] - 1) * (2s[2] - 1))
    2 ≤ move ≤ 3 && (ΔE += 2 * (2s[2] - 1) * (2s[3] - 1))
    (move == 1 || move == 3) && (ΔE += 2 * (2s[3] - 1) * (2s[1] - 1))

    return ΔE
end

#neighbors(X::GraphThreeSpin, i::Int) = return i == 1 ? (2,) : i == 2 ? (1, 3) : (2,)
neighbors(X::GraphThreeSpin, i::Int) = return (mod1(i-1, 3), mod1(i+1,3))
allΔE(::Type{GraphThreeSpin}) = (0, 4)




end
