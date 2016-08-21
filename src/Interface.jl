module Interface

using ExtractMacro
using ..Common

export Config, AbstractGraph, SimpleGraph, DiscrGraph, DoubleGraph, spinflip!,
       energy, delta_energy, neighbors, getN, allΔE, discr_graph,
       delta_energy_residual, update_cache!, update_cache_residual!

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

#spinflip!(C::Config, move::Int) = (C.s[move] $= 1)
spinflip!(C::Config, move::Int) = unsafe_bitflip!(C.s, move)

"""
    AbstractGraph{ET<:Real}

An abstract type representing an Ising spin model. The `ET` parameter
is the type returned by the [`energy`](@ref) and [`delta_energy`](@ref)
functions.

See also [`DiscrGraph`](@ref) and [`DoubleGraph`](@ref).
"""
abstract AbstractGraph{ET<:Real}

"""
    update_cache!(X::AbstractGraph, C::Config, move::Int)

A function which is called every time a spin is flipped. This may happen:

1. when a move is accepted, in [`standardMC`](@ref) and [`rrrMC`](@ref)
2. when a move is attempted to evaluate the effect on the neighbors, in [`rrrMC`](@ref) and [`bklMC`](@ref).

`move` is the spin index. By default, this function does nothing, but it may be overloaded
by particular graph types.

When `X` is a [`DoubleGraph`](@ref), there is a default implementation which first calls
`update_cache!` on [`discr_graph`](@ref)`(X)`, then
calls [`update_cache_residual!`](@ref) on `X`.

*Note*: this function is always invoked *after* the flip has been performed, unlike in [`delta_energy`](@ref)
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
beginning of [`standardMC`](@ref), [`rrrMC`](@ref) and [`bklMC`](@ref). Subsequently,
[`delta_energy`](@ref) is used instead.

All graphs must implement this function.

It should also be used to initialize/reset the cache for a given graph, if any (see [`update_cache!`](@ref)).
"""
energy(::AbstractGraph, ::Config) = error("not implemented")

"""
    delta_energy(X::AbstractGraph, C::Config, move::Int)

Returns the energy difference that would be associated to flipping the spin `move`.

A default fallback implementation based on `energy` is provided, to be used for debugging,
but having an efficient implementation for each graph is critical for performance.

*Note*: when `X` is a [`DiscrGraph`](@ref), the absolute value of the result must be contained in the
tuple returned by [`allΔE`](@ref) – no approximations are allowed, and missing values will cause crashes
(unless Julia is run with the `--check-bounds=yes` option, in which case they will cause errors).

*Note*: this function is always invoked *before* performing the flip, unlike in [`update_cache!`](@ref)
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

Returns the number of spins for a graph. Defaults to just `X.N`.
"""
getN(X::AbstractGraph) = X.N


"""
    SimpleGraph{ET} <: AbstractGraph{ET}

An abstract type representing a generic graph. This can only be used with [`standardMC`](@ref),
not with [`rrrMC`](@ref) or [`bklMC`](@ref).

The `ET` parameter is the type returned by [`energy`](@ref) and [`delta_energy`](@ref).
"""
abstract SimpleGraph{ET} <: AbstractGraph{ET}


"""
    DiscrGraph{ET} <: AbstractGraph{ET}

An abstract type representing a graph in which the [`delta_energy`](@ref) values
produced when flipping a spin belong to a finite discrete set, and thus can be
sampled efficiently with [`rrrMC`](@ref) or [`bklMC`](@ref).

The `ET` parameter is the type returned by [`energy`](@ref) and [`delta_energy`](@ref).

It is also used internally in [`DoubleGraph`](@ref).

See also [`neighbors`](@ref) and [`allΔE`](@ref).
"""
abstract DiscrGraph{ET} <: AbstractGraph{ET}

"""
    neighbors(X::DiscrGraph, i::Int)

Returns an iterable with all the neighbors of spin `i`. This is required
by [`rrrMC`](@ref) and [`bklMC`](@ref) since those methods need to evaluate
the effect of flipping a spin on its neighbors' delta-energy classes.

For performance reasons, it is best if the returned value is stack-allocated
rather than heap-allocated, e.g. it is better to return a `Tuple` than a `Vector`.
"""
neighbors(::DiscrGraph, i::Int) = error("not implemented")

"""
    allΔE{P<:DiscrGraph}(::Type{P})

Returns a tuple of all possible *positive* values that can be returned
by [`delta_energy`](@ref). This must be implemented by all `DiscrGraph`
objects in order to use [`rrrMC`](@ref) or [`bklMC`](@ref).

For performance reasons, it is best if the result can be computed
from the type of the graph alone (possibly using a generated function).
"""
allΔE{P<:DiscrGraph}(::Type{P}) = error("not implemented")
allΔE(X::DiscrGraph) = allΔE(typeof(X))

"""
    DoubleGraph{ET} <: AbstractGraph{ET}

An abstract type representing a graph in which the energy is the sum of two
contributions, one of which can be encoded in a [`DiscrGraph`](@ref) type.
This allows [`rrrMC`](@ref) to sample values more efficiently.

The `ET` parameter is the type returned by the [`energy`](@ref) and
[`delta_energy`](@ref) functions. Note that it can be different from
the type of the internal `DiscrGraph` object (e.g., one can have
a `DiscrGraph{Int}` object inside a `DoubleGraph{Float64}` object).

See also [`discr_graph`](@ref), [`delta_energy_residual`](@ref) and
[`update_cache_residual!`](@ref).
"""
abstract DoubleGraph{ET} <: AbstractGraph{ET}

"""
    discr_graph(X::DoubleGraph)

Returns the internal [`DiscrGraph`](@ref) used by the given
[`DoubleGraph`](@ref). The default is simply `X.X0`.
"""
discr_graph(X::DoubleGraph) = X.X0
discr_graph(X::DiscrGraph) = X

"""
    delta_energy_residual(X::DoubleGraph, C::Config, move::Int)

Returns the residual part of the energy difference produced if the spin `move` would
be flipped, excluding the contribution from the internal [`DiscrGraph`](@ref)
(see [`discr_graph`](@ref)).

See also [`delta_energy`](@ref). There is a default fallback implementation, but
it should be overloaded for efficiency.
"""

delta_energy_residual(X::DoubleGraph, C::Config, move::Int) =
    delta_energy(X, C, move) - delta_energy(discr_graph(X), C, move)

"""
    update_cache_residual!(X::DoubleGraph, C::Config, move::Int)

Called internally by the default [`update_cache!`](@ref) when the argument is a [`DoubleGraph`](@ref).
Can be useful to overload this if the residual part of the graph has an indipendent cache.

By default, it does nothing.
"""
update_cache_residual!(X::DoubleGraph, C::Config, move::Int) = nothing

function update_cache!(X::DoubleGraph, C::Config, move::Int)
    update_cache!(discr_graph(X), C, move)
    update_cache_residual!(X, C, move)
end

allΔE(X::DoubleGraph) = allΔE(typeof(discr_graph(X)))

end # module
