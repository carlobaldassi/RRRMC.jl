# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module Interface

using Compat, ExtractMacro
using ..Common

if isdefined(Main, :Documenter)
# this is silly but it's required for correct cross-linking in docstrings, apparently
using ...RRRMC
end

export Config, AbstractGraph, SimpleGraph, DiscrGraph, SingleGraph, DoubleGraph,
       spinflip!, energy, delta_energy, neighbors, getN, allΔE, inner_graph,
       delta_energy_residual, update_cache!, update_cache_residual!,
       cenergy, distances

import Base: length, copy!, copy, ==

immutable Config
    N::Int
    s::BitVector
    function Config(N::Integer; init::Bool = true)
        s = BitArray(N)
        init && rand!(s)
        return new(N, s)
    end
end

@doc """
    Config(N::Integer)

The object storing the configuration for an Ising model. Although
the spin values are \$σ_i ∈ {-1,1}\$, internally they are stored in
a `BitArray`, in the type field `s`, so that to obtain the real
value one needs to perform the transformation \$σ_i = 2s_i - 1\$.
""" -> Config(N::Integer)

length(C::Config) = C.N

#spinflip!(C::Config, move::Int) = (C.s[move] $= 1)
spinflip!(C::Config, move::Int) = unsafe_bitflip!(C.s, move)

function copy!(dst::Config, src::Config)
    length(dst) == length(src) || throw(ArgumentError("Incompatible Config lengths, dst=$(length(dst)) src=$(length(src))"))
    copy!(dst.s, src.s)
    return dst
end

copy(C::Config) = copy!(Config(length(C), init=false), C)

(==)(C1::Config, C2::Config) = length(C1) == length(C2) && C1.s == C2.s

"""
    AbstractGraph{ET<:Real}

An abstract type representing an Ising spin model. The `ET` parameter
is the type returned by the [`energy`](@ref) and [`delta_energy`](@ref)
functions.

See also [`SimpleGraph`](@ref), [`DiscrGraph`](@ref), [`SingleGraph`](@ref)
and [`DoubleGraph`](@ref).
"""
@compat abstract type AbstractGraph{ET<:Real} end

"""
    update_cache!(X::AbstractGraph, C::Config, move::Int)

A function which is called every time a spin is flipped. This may happen:

1. when a move is accepted, in [`standardMC`](@ref), [`rrrMC`](@ref), [`bklMC`](@ref) and
   [`wtmMC`](@ref).
2. when a move is attempted to evaluate the effect on the neighbors, in [`rrrMC`](@ref).

`move` is the spin index. By default, this function does nothing, but it may be overloaded
by particular graph types.

When `X` is a [`DoubleGraph`](@ref), there is a default implementation which first calls
`update_cache!` on [`inner_graph`](@ref)`(X)`, then
calls [`update_cache_residual!`](@ref) on `X`.

*Note:* this function is always invoked *after* the flip has been performed, unlike in [`delta_energy`](@ref)
and [`delta_energy_residual`](@ref).
"""
update_cache!(X::AbstractGraph, C::Config, move::Int) = nothing

function spinflip!(X::AbstractGraph, C::Config, move::Int)
    spinflip!(C, move)
    update_cache!(X, C, move)
end

"""
    energy(X::AbstractGraph, C::Config)

Returns the energy of graph `X` in the configuration `C`. This is always invoked at the
beginning of [`standardMC`](@ref), [`rrrMC`](@ref), [`bklMC`](@ref) and [`wtmMC`](@ref).
Subsequently, [`delta_energy`](@ref) is used instead.

All graphs must implement this function.

It *must* also be used to initialize/reset the cache for a given graph, if any (see [`update_cache!`](@ref)).
"""
energy(::AbstractGraph, ::Config) = error("not implemented")

"""
    delta_energy(X::AbstractGraph, C::Config, move::Int)

Returns the energy difference that would be associated to flipping the spin `move`.

A default fallback implementation based on `energy` is provided, to be used for debugging,
but having an efficient implementation for each graph is critical for performance.

*Note:* when `X` is a [`DiscrGraph`](@ref), the absolute value of the result must be contained in the
tuple returned by [`allΔE`](@ref) – no approximations are allowed, and missing values will cause crashes
(unless Julia is run with the `--check-bounds=yes` option, in which case they will cause errors).

*Note:* this function is always invoked *before* performing the flip, unlike in [`update_cache!`](@ref)
and [`update_cache_residual!`](@ref).
"""
function delta_energy(X::AbstractGraph, C::Config, move::Int)
    @extract C : s
    oldn = energy(X, C)
    s[move] $= 1
    newn = energy(X, C)
    s[move] $= 1
    Δn0 = newn - oldn
    return Δn0
end

"""
    getN(X::AbstractGraph)

Returns the number of spins for a graph. The default implementation just returns `X.N`.
"""
getN(X::AbstractGraph) = X.N

"""
    neighbors(X::AbstractGraph, i::Int)

Returns an iterable with all the neighbors of spin `i`. This is required
by [`rrrMC`](@ref), [`bklMC`](@ref) and [`wtmMC`](@ref) since those methods need to
evaluate the effect of flipping a spin on its neighbors' local fields. It is
not required by [`standardMC`](@ref).

For performance reasons, it is best if the returned value is stack-allocated
rather than heap-allocated, e.g. it is better to return an `NTuple` than a `Vector`.
"""
neighbors(::AbstractGraph, i::Int) = error("not implemented")



"""
    SimpleGraph{ET} <: AbstractGraph{ET}

An abstract type representing a generic graph.

The `ET` parameter is the type returned by [`energy`](@ref) and [`delta_energy`](@ref).
"""
@compat abstract type SimpleGraph{ET} <: AbstractGraph{ET} end


"""
    DiscrGraph{ET} <: AbstractGraph{ET}

An abstract type representing a graph in which the [`delta_energy`](@ref) values
produced when flipping a spin belong to a finite discrete set, and thus can be
sampled more efficiently with [`rrrMC`](@ref) or [`bklMC`](@ref) if the set is
small.

The `ET` parameter is the type returned by [`energy`](@ref) and [`delta_energy`](@ref).

See also [`allΔE`](@ref).
"""
@compat abstract type DiscrGraph{ET} <: AbstractGraph{ET} end

"""
    allΔE{P<:DiscrGraph}(::Type{P})
    allΔE(X::DiscrGraph)

Returns a tuple of all possible *non-negative* values that can be returned
by [`delta_energy`](@ref). This must be implemented by all `DiscrGraph`
objects in order to use [`rrrMC`](@ref) or [`bklMC`](@ref).

The second form above, called directly on the `DiscrGraph` object, takes
precedence; however, for performance reasons, it is best if the result can
be computed from the type of the graph alone (possibly using a generated
function), i.e. it is best to implement the first definition above — the
second one then falls back automatically to the first.
"""
allΔE{P<:DiscrGraph}(::Type{P}) = error("not implemented")
allΔE(X::DiscrGraph) = allΔE(typeof(X))

"""
    SingleGraph{ET}

A type alias representing either a [`SimpleGraph`](@ref) or a
[`DiscrGraph{ET}`](@ref). See also [`DoubleGraph`](@ref).
"""
@compat const SingleGraph{ET} = Union{SimpleGraph{ET},DiscrGraph{ET}}

"""
    DoubleGraph{GT<:SingleGraph,ET} <: AbstractGraph{ET}

An abstract type representing a graph in which the energy is the sum of two
contributions, one of which is encoded in a graph of type `GT` (see
[`SingleGraph`](@ref)). This allows [`rrrMC`](@ref) to sample values
more efficiently.

The `ET` parameter is the type returned by the [`energy`](@ref) and
[`delta_energy`](@ref) functions. Note that it can be different from
the energy type of the internal `GT` object (e.g., one can have
a `DoubleGraph{DiscrGraph{Int},Float64}` object).

*Note:* When you declare a type as subtype of this, `GT` should *not* be the
concrete type of the inner graph, but either `SimpleGraph{T}` or
`DiscrGraph{T}` for some `T`.

See also [`inner_graph`](@ref), [`delta_energy_residual`](@ref) and
[`update_cache_residual!`](@ref).
"""
@compat abstract type DoubleGraph{GT,ET} <: AbstractGraph{ET} end

"""
    inner_graph(X::DoubleGraph)

Returns the internal graph used by the given [`DoubleGraph`](@ref).
The default implementation simply returns `X.X0`.
"""
inner_graph(X::DoubleGraph) = X.X0
inner_graph(X::SingleGraph) = X

"""
    delta_energy_residual(X::DoubleGraph, C::Config, move::Int)

Returns the residual part of the energy difference produced if the spin `move` would
be flipped, excluding the contribution from the internal [`SimpleGraph`](@ref)
(see [`inner_graph`](@ref)).

See also [`delta_energy`](@ref). There is a default fallback implementation, but
it should be overloaded for efficiency.
"""

delta_energy_residual(X::DoubleGraph, C::Config, move::Int) =
    delta_energy(X, C, move) - delta_energy(inner_graph(X), C, move)

"""
    update_cache_residual!(X::DoubleGraph, C::Config, move::Int)

Called internally by the default [`update_cache!`](@ref) when the argument is a [`DoubleGraph`](@ref).
Can be useful to overload this if the residual part of the graph has an indipendent cache.

By default, it does nothing.
"""
update_cache_residual!(X::DoubleGraph, C::Config, move::Int) = nothing

function update_cache!(X::DoubleGraph, C::Config, move::Int)
    update_cache!(inner_graph(X), C, move)
    update_cache_residual!(X, C, move)
end

allΔE{ET,GT<:DiscrGraph}(X::DoubleGraph{ET,GT}) = allΔE(typeof(inner_graph(X)))

"""
    cenergy(X::DoubleGraph)

Returns the individual energy (as defined by the original model)
of the reference configuration in a [`GraphLocalEntropy`](@ref) graph or a
[`GraphTopologicalLocalEntropy`](@ref) graph.
"""
cenergy(::AbstractGraph) = error("not implemented")

"""
    distances(X::DoubleGraph)

Returns the matrix of the Hamming distances between replicas, in graphs with
replicas: [`GraphQuant`](@ref), [`GraphLocalEntropy`](@ref), [`GraphTopologicalLocalEntropy`](@ref)
or [`GraphRobustEnsemble`](@ref).
"""
distances(::AbstractGraph) = error("not implemented")

end # module
