# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module Empty

using ExtractMacro
using ..Interface
using ..Common
using ..QT

export GraphEmpty

import ..Interface: energy, delta_energy, neighbors

struct GraphEmpty <: SimpleGraph{Int}
    N::Int
    GraphEmpty(N::Integer) = return new(N)
end

@doc """
    GraphEmpty(N::Integer) <: SimpleGraph{Int}

A trivial `SimpleGraph` type with `N` free, non-interacting spins.
(The energy is always `0`.)

Only useful for testing/debugging purposes.
""" -> GraphEmpty()

energy(X::GraphEmpty, C::Config) = 0
delta_energy(X::GraphEmpty, C::Config, move::Int) = 0

neighbors(X::GraphEmpty, i::Int) = ()

end # module
