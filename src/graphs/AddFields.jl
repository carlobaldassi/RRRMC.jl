# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module AddFields

using ExtractMacro
using ..Interface
using ..Common

using ...RRRMC # this is silly but it's required for correct cross-linking in docstrings, apparently

export GraphAF, GraphAddFields, GraphAddSubFields

import ..Interface: energy, delta_energy, neighbors,
                    update_cache!, delta_energy_residual

struct GraphAF{ET} <: SimpleGraph{ET}
    N::Int
    fields::Vector{ET}
    function GraphAF(fields::Vector{ET}) where {ET}
        N = length(fields)
        return new{ET}(N, fields)
    end
end

@doc """
    GraphAF(fields::Vector) <: SimpleGraph

An auxiliary `SimpleGraph` used to implement additional fields on variables.

It is only useful when implementing other graph types; see [`GraphAddFields`](@ref).
""" -> GraphAF(fields::Vector)

function energy(X::GraphAF{ET}, C::Config) where {ET}
    @assert X.N == C.N
    @extract X : N fields
    @extract C : s
    E = zero(ET)
    for i = 1:N
        σi = 2s[i] - 1
        E += fields[i] * σi
    end
    return E
end

function delta_energy(X::GraphAF, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : σ = 2s[move]-1
    @extract X : h = fields[move]

    return -2 * σ * h
end

neighbors(X::GraphAF, i::Int) = ()

# Add external fields to any AbstractGraph

mutable struct GraphAddFields{ET,G<:AbstractGraph} <: DoubleGraph{SimpleGraph{ET},ET}
    N::Int
    X0::GraphAF{ET}
    X1::G
    function GraphAddFields(fields::Vector{ET}, g::G) where {ET,G}
        X0 = GraphAF(fields)
        N = X0.N
        getN(g) == N || throw(ArgumentError("incompatible length, fields size=$N graph size=$(getN(g))"))
        return new{ET,G}(N, X0, g)
    end
end

@doc """
    GraphAddFields(fields::Vector, g::AbstractGraph) <: DoubleGraph

A `DoubleGraph` that implements a model in which additional `fields` are added to the
graph `g`. The additional fields can then be taken care of efficiently with [`rrrMC`](@ref).
""" -> GraphAddFields(fields::Vector, g::AbstractGraph)

function update_cache!(X::GraphAddFields, C::Config, move::Int)
    update_cache!(X.X1, C, move)
end

energy(X::GraphAddFields, C::Config) = energy(X.X0, C) + energy(X.X1, C)

delta_energy_residual(X::GraphAddFields, C::Config, move::Int) = delta_energy(X.X1, C, move)

function delta_energy(X::GraphAddFields, C::Config, move::Int)
    return delta_energy(X.X0, C, move) +
           delta_energy_residual(X, C, move)
end

neighbors(X::GraphAddFields, i::Int) = neighbors(X.X1, i)



mutable struct GraphAddSubFields{ET,G<:AbstractGraph} <: DoubleGraph{SimpleGraph{ET},ET}
    N::Int
    X0::GraphAF{ET}
    X1::G
    function GraphAddSubFields(fields::Vector{ET}, g::G) where {ET,G}
        X0 = GraphAF(fields)
        N = X0.N
        getN(g) == N || throw(ArgumentError("incompatible length, fields size=$N graph size=$(getN(g))"))
        return new{ET,G}(N, X0, g)
    end
end

@doc """
    GraphAddSubFields(fields::Vector, g::AbstractGraph) <: DoubleGraph

A `DoubleGraph` that implements a model in which the given `fields` are added and subtracted to the
graph `g`. The fields can then be taken care of efficiently with [`rrrMC`](@ref).
""" -> GraphAddSubFields(fields::Vector, g::AbstractGraph)

function update_cache!(X::GraphAddSubFields, C::Config, move::Int)
    update_cache!(X.X1, C, move)
end

energy(X::GraphAddSubFields, C::Config) = energy(X.X1, C)

delta_energy_residual(X::GraphAddSubFields, C::Config, move::Int) = delta_energy(X.X1, C, move) - delta_energy(X.X0, C, move)

delta_energy(X::GraphAddSubFields, C::Config, move::Int) = delta_energy(X.X1, C, move)

neighbors(X::GraphAddSubFields, i::Int) = neighbors(X.X1, i)

end # module
