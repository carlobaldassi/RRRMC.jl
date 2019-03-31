# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module RE

using ExtractMacro
using ..Interface
using ..Common

using ...RRRMC # this is silly but it's required for correct cross-linking in docstrings, apparently

export GraphRE, GraphRobustEnsemble, REenergies

import ..Interface: energy, delta_energy, neighbors, allΔE,
                    update_cache!, delta_energy_residual

import Base: iterate, length, eltype

function logcoshratio(a, b)
    #return log(cosh(a) / cosh(b))
    #log((exp(a) + exp(-a)) / (exp(b) + exp(-b)))
    a = abs(a)
    b = abs(b)
    return a - b + (log1p(exp(-2a)) - log1p(exp(-2b)))
end

fk(μ̄::Integer, γ::Real, β::Real) = logcoshratio(γ * (μ̄ + 1), γ * (μ̄ - 1)) / β

mutable struct GraphRE{M,γ,β} <: DiscrGraph{Float64}
    N::Int
    Nk::Int
    μ::IVec
    cache::LocalFields{Float64}

    function GraphRE{M,γ,β}(N::Integer, μ::Union{Nothing,IVec} = nothing) where {M,γ,β}
        isa(M, Int) || throw(ArgumentError("invalid parameter M, expected Int, given: $(typeof(M))"))
        M > 2 || throw(ArgumentError("M must be greater than 2, given: $M"))
        isa(β, Float64) || throw(ArgumentError("invalid parameter β, expected Float64, given: $(typeof(β))"))
        isa(γ, Float64) || throw(ArgumentError("invalid parameter γ, expected Float64, given: $(typeof(γ))"))
        N % M == 0 || throw(ArgumentError("N must be divisible by M, given: N=$N M=$M"))
        Nk = N ÷ M

        μ::IVec = (μ ≡ nothing ? zeros(Int, Nk) : μ)
        @assert length(μ) == Nk
        cache = LocalFields{Float64}(N)
        return new{M,γ,β}(N, Nk, μ, cache)
    end
end

@doc """
    GraphRE{M,γ,β}(N::Integer) <: DiscrGraph

An auxiliary `DiscrGraph` used to implement the interactions in the
Robust Ensemble.

It is only useful when implementing other graph types; see [`GraphRobustEnsemble`](@ref).
""" -> GraphRE{M,γ,β}(N::Integer)

GraphRE(X::GraphRE{M,oldγ,β}, newγ::Float64) where {M,oldγ,β} = GraphRE{M,newγ,β}(X.N, X.μ)

@generated function ΔElist(::Type{GraphRE{M,γ,β}}) where {M,γ,β}
    Expr(:tuple, ntuple(d->fk(2*(d - 1 - (M-1) >>> 0x1) - iseven(M), γ, β), M)...)
end
lstind(μ̄::Int, M::Int) = (μ̄ + M-1) >>> 0x1 + 1

function getk(X::GraphRE{M,γ,β}, μ̄::Int) where {M,γ,β}
    @inbounds k = ΔElist(GraphRE{M,γ,β})[lstind(μ̄, M)]
    return k
end

function energy(X::GraphRE{M,γ,β}, C::Config) where {M,γ,β}
    # @assert X.N == C.N
    @extract X : Nk μ cache
    @extract cache : lfields lfields_last
    @extract C : s

    fill!(μ, 0)
    j = 0
    for i = 1:Nk, k = 1:M
        j += 1
        # @assert j == k + (i-1) * M
        σj = 2s[j] - 1
        μ[i] += σj
    end

    n = 0.0
    for i = 1:Nk
        n -= log(2cosh(γ * μ[i])) / β
    end

    j = 0
    for i = 1:Nk, k = 1:M
        j += 1
        # @assert j == k + (i-1) * M
        σj = 2s[j] - 1
        μ̄ = μ[i] - σj
        # kj = fk(μ̄, γ, β)
        kj = getk(X, μ̄)
        # @assert kj == fk(μ̄, γ, β)
        lfields[j] = σj * kj
    end
    cache.move_last = 0
    fill!(lfields_last, 0.0)
    return n
end

function kinterval(move::Integer, M::Integer)
    j0 = move - ((move-1) % M)
    j1 = j0 + M - 1
    return j0:j1
end

function update_cache!(X::GraphRE{M,γ,β}, C::Config, move::Int) where {M,γ,β}
    # @assert X.N == C.N
    # @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : Nk μ cache

    @inbounds σx = 2s[move] - 1
    i = (move-1) ÷ M + 1

    Ux = kinterval(move, M)

    @extract cache : lfields lfields_last move_last
    if move_last == move
        @inbounds begin
            #for y in neighbors(X, move)
            for y in Ux
                lfields[y], lfields_last[y] = lfields_last[y], lfields[y]
            end
            #lfields[move] = -lfields[move]
            #lfields_last[move] = -lfields_last[move]
        end
        μ[i] += 2σx
        return
    end

    @inbounds begin
        μnew = μ[i] + 2σx
        μ[i] = μnew
        #for y in neighbors(X, move)
            # @assert y ≠ move
        for y in Ux
            σy = 2s[y] - 1
            μ̄ = μnew - σy
            # @assert -M+1 ≤ μ̄ ≤ M-1    μ̄
            ky = getk(X, μ̄)
            # @assert ky == fk(μ̄, γ, β)
            lfields_last[y] = lfields[y]
            lfields[y] = σy * ky
        end
        #lfm = lfields[move]
        #lfields_last[move] = lfm
        #lfields[move] = -lfm
    end
    cache.move_last = move

    # lfields_bk = copy(lfields)
    # energy(X, C)
    # lfields_bk ≠ lfields && @show move hcat(lfields,lfields_bk)
    # @assert lfields_bk == lfields

    return
end

@inline function delta_energy(X::GraphRE, C::Config, move::Int)
    # @assert X.N == C.N
    # @assert 1 ≤ move ≤ C.N
    @extract X : cache
    @extract cache : lfields

    @inbounds Δ = lfields[move]
    return Δ
end

struct CavityRange
    j0::Int
    j1::Int
    jX::Int
    function CavityRange(j0::Integer, j1::Integer, jX::Integer)
        j0 ≤ jX ≤ j1 || throw(ArgumentError("invalid CavityRange parameters, expected j0≤jX≤j1, given: j0=$j0, j1=$j1, jX=$X"))
        return new(j0, j1, jX)
    end
end

function iterate(crange::CavityRange, j = crange.j0 + (crange.jX == crange.j0))
    j > crange.j1 && return nothing
    return j, j + 1 + (j == crange.jX - 1)
end

# start(crange::CavityRange) = crange.j0 + (crange.jX == crange.j0)
# done(crange::CavityRange, j) = j > crange.j1
# @inline function next(crange::CavityRange, j)
#     @extract crange : j0 j1 jX
#     # @assert j ≠ jX
#     nj = j + 1
#     nj += (nj == jX)
#     return (j, nj)
# end
length(crange::CavityRange) = crange.j1 - crange.j0
eltype(::Type{CavityRange}) = Int

@inline function neighbors(X::GraphRE{M}, j::Int) where {M}
    j0 = j - ((j-1) % M)
    j1 = j0 + M - 1
    return CavityRange(j0, j1, j)
end

@generated function allΔE(::Type{GraphRE{M,γ,β}}) where {M,γ,β}
    K = M - 1
    iseven(K) ? Expr(:tuple, ntuple(d->fk(2*(d-1), γ, β), K÷2+1)...) :
                Expr(:tuple, ntuple(d->fk(2d-1, γ, β), (K+1)÷2)...)
end

# Replicate an existsing graph

mutable struct GraphRobustEnsemble{M,γ,β,G<:AbstractGraph} <: DoubleGraph{DiscrGraph{Float64},Float64}
    N::Int
    Nk::Int
    X0::GraphRE{M,γ,β}
    X1::Vector{G}
    C1::Vector{Config}
    function GraphRobustEnsemble{M,γ,β,G}(N::Integer, g0::G, Gconstr, args...) where {M,γ,β,G}
        X0 = GraphRE{M,γ,β}(N)
        Nk = X0.Nk
        X1 = Array{G}(undef, M)
        X1[1] = g0
        for k = 2:M
            X1[k] = Gconstr(args...)
        end
        C1 = [Config(Nk, init=false) for k = 1:M]
        return new{M,γ,β,G}(N, Nk, X0, X1, C1)
    end
end

"""
    GraphRobustEnsemble(N::Integer, M::Integer, γ::Float64, β::Float64, Gconstr, args...) <: DoubleGraph

A `DoubleGraph` which implements a "Robust Ensemble" model, given any other Ising model previously defined.
This kind of graph can be simulated efficiently with [`rrrMC`](@ref).

`N` is the number of spins, `M` the number of replicas, `γ` the interaction strength,
`β` the inverse temperature. `GConstr` is the (original) graph constructor, and `args` the
arguments to the contructor.

This is similar to [`GraphLocalEntropy`](@ref RRRMC.GraphLocalEntropy), but the reference configuration is traced out.

See also [`REenergies`](@ref).
"""
function GraphRobustEnsemble(Nk::Integer, M::Integer, γ::Float64, β::Float64, Gconstr, args...)
    g0 = Gconstr(args...)
    G = typeof(g0)
    return GraphRobustEnsemble{M,γ,β,G}(Nk * M, g0, Gconstr, args...)
end

function update_cache!(X::GraphRobustEnsemble{M}, C::Config, move::Int) where {M}
    @extract X : X0 X1 C1
    k = mod1(move, M)
    i = (move - 1) ÷ M + 1

    spinflip!(X1[k], C1[k], i)

    update_cache!(X0, C, move)
end

function energy(X::GraphRobustEnsemble{M}, C::Config) where {M}
    # @assert X.N == C.N
    @extract X : Nk X0 X1 C1
    @extract C : s

    E = energy(X0, C)

    for k = 1:M
        s1 = C1[k].s
        for (i,j) = enumerate(k:M:(k + M * (Nk-1)))
            s1[i] = s[j]
        end
        E += energy(X1[k], C1[k])
    end

    return E
end

"""
    REenergies(X::GraphRobustEnsemble)

Returns a Vector with the individual energy (as defined by the original model)
of each replica in a [`GraphRobustEnsemble`](@ref) graph.
"""
function REenergies(X::GraphRobustEnsemble{M}) where {M}
    @extract X : X1 C1

    Es = zeros(M)

    for k = 1:M
        Es[k] = energy(X1[k], C1[k])
    end

    return Es
end

function delta_energy_residual(X::GraphRobustEnsemble{M}, C::Config, move::Int) where {M}
    @extract X : X1 C1

    k = mod1(move, M)
    i = (move - 1) ÷ M + 1

    return delta_energy(X1[k], C1[k], i)
end

function delta_energy(X::GraphRobustEnsemble, C::Config, move::Int)
    return delta_energy(X.X0, C, move) +
           delta_energy_residual(X, C, move)
end

# This is inefficient. On the other hand, it is
# basically only written for testing purposes...
function neighbors(X::GraphRobustEnsemble{M}, i::Int) where {M}
    @extract X : X0 X1
    jts = neighbors(X0, i)

    k = mod1(i, M)
    i1 = (i - 1) ÷ M + 1

    is = neighbors(X1[k], i1)

    return tuple(jts..., map(i->((i-1)*M+k), is)...)
end

end # module
