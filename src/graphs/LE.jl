# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module LE

using ExtractMacro
using ..Interface
using ..Common

if isdefined(Main, :Documenter)
using ...RRRMC # this is silly but it's required for correct cross-linking in docstrings, apparently
end

export GraphLE, GraphLocalEntropy, LEenergies

import ..Interface: energy, delta_energy, neighbors, allΔE,
                    update_cache!, delta_energy_residual,
                    cenergy, distances

type GraphLE{M,γT} <: DiscrGraph{Float64}
    N::Int
    Nk::Int
    cache::LocalFields{Int}

    @inner {M,γT} function GraphLE(N::Integer)
        isa(M, Int) || throw(ArgumentError("invalid parameter M, expected Int, given: $(typeof(M))"))
        M > 2 || throw(ArgumentError("M must be greater than 2, given: $M"))
        isa(γT, Float64) || throw(ArgumentError("invalid parameter γT, expected Float64, given: $(typeof(γT))"))
        N % (M + 1) == 0 || throw(ArgumentError("N must be divisible by M+1, given: N=$N M=$M"))
        Nk = N ÷ (M+1)

        cache = LocalFields{Int}(N)
        return new(N, Nk, cache)
    end
end

@doc """
    GraphLE{M,γT}(N::Integer) <: DiscrGraph

An auxiliary `DiscrGraph` used to implement the interactions in the
Local Entropy Ensemble.

It is only useful when implementing other graph types; see [`GraphLocalEntropy`](@ref).
""" -> GraphLE{M,γT}(N::Integer)

GraphLE{M,oldγ}(X::GraphLE{M,oldγ}, newγ::Float64) = GraphLE{M,newγ}(X.N)

# @generated function ΔElist{M,γT}(::Type{GraphLE{M,γT}})
#     Expr(:tuple, ntuple(d->fk(2*(d - 1 - (M-1) >>> 0x1) - iseven(M), γT), M)...)
# end
# lstind(μ̄::Int, M::Int) = (μ̄ + M-1) >>> 0x1 + 1
#
# function getk{M,γT}(X::GraphLE{M,γT}, μ̄::Int)
#     @inbounds k = ΔElist(GraphLE{M,γT})[lstind(μ̄, M)]
#     return k
# end

function energy{M,γT}(X::GraphLE{M,γT}, C::Config)
    # @assert X.N == C.N
    @extract X : Nk cache
    @extract cache : lfields lfields_last
    @extract C : s

    n = 0
    j = 0
    for i = 1:Nk
        j += 1
        # @assert j == 1 + (i-1) * (M+1)
        jc = j
        σc = 2s[jc] - 1

        μ = 0
        for k = 2:(M+1)
            j += 1
            # @assert j == k + (i-1) * (M+1)
            σj = 2s[j] - 1
            lfields[j] = σc * σj
            μ += σj
        end
        f = σc * μ
        lfields[jc] = f
        n -= f
    end
    cache.move_last = 0
    fill!(lfields_last, 0.0)
    return n * γT
end

function kinterval(move::Integer, M::Integer)
    j0 = move - ((move-1) % (M+1))
    j1 = j0 + M
    return j0:j1
end

function update_cache!{M,γT}(X::GraphLE{M,γT}, C::Config, move::Int)
    # @assert X.N == C.N
    # @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : Nk cache

    k = mod1(move, M+1)
    if k ≠ 1
        @inbounds σx = 2s[move] - 1
        i = (move-1) ÷ (M+1) + 1
        jc = 1 + (i-1) * (M+1)
    else
        Ux = kinterval(move, M)
    end

    @extract cache : lfields lfields_last move_last
    if move_last == move
        @inbounds begin
            #for y in neighbors(X, move)
            if k ≠ 1
                lfields[jc], lfields_last[jc] = lfields_last[jc], lfields[jc]
                lfields[move], lfields_last[move] = lfields_last[move], lfields[move]
            else
                for y in Ux
                    lfields[y], lfields_last[y] = lfields_last[y], lfields[y]
                end
            end
            #lfields[move] = -lfields[move]
            #lfields_last[move] = -lfields_last[move]
        end
        return
    end

    @inbounds begin
        if k ≠ 1
            #jc = 1 + (i-1) * (M+1)
            σc = 2s[jc] - 1
            lfc = lfields[jc]
            lfm = lfields[move]
            lfields_last[jc] = lfc
            lfields_last[move] = lfm
            lfields[jc] = lfc + 2 * (σc * σx)
            lfields[move] = -lfm
        else
            for y in Ux
                lfy = lfields[y]
                lfields_last[y] = lfy
                lfields[y] = -lfy
            end
        end
        #lfm = lfields[move]
        #lfields_last[move] = lfm
        #lfields[move] = -lfm
    end
    cache.move_last = move

    # lfields_bk = copy(lfields)
    # energy(X, C)
    # lfields_bk ≠ lfields && @show move hcat(lfields,lfields_bk) find(lfields .≠ lfields_bk)
    # @assert lfields_bk == lfields

    return
end

@inline function delta_energy{M,γT}(X::GraphLE{M,γT}, C::Config, move::Int)
    # @assert X.N == C.N
    # @assert 1 ≤ move ≤ C.N
    @extract X : cache
    @extract cache : lfields

    @inbounds Δ = 2γT * lfields[move]
    return Δ
end

@inline function neighbors{M}(X::GraphLE{M}, j::Int)
    r = (j-1) % (M+1)
    if r == 0
        return (j+1):(j+M)
    else
        j0 = j - r
        return j0:j0
    end
end

@generated function allΔE{M,γT}(::Type{GraphLE{M,γT}})
    iseven(M) ? Expr(:tuple, insert!([4*d*abs(γT) for d = 0:(M÷2)], 2, 2abs(γT))...) :
                Expr(:tuple, [2*(2d-1)*abs(γT) for d = 1:((M+1)÷2)]...)
end

# Replicate an existsing graph

type GraphLocalEntropy{M,γT,G<:AbstractGraph} <: DoubleGraph{DiscrGraph{Float64},Float64}
    N::Int
    Nk::Int
    X0::GraphLE{M,γT}
    Xc::G
    X1::Vector{G}
    Cc::Config
    C1::Vector{Config}
    @inner {M,γT,G} function GraphLocalEntropy(N::Integer, g0::G, Gconstr, args...)
        X0 = GraphLE{M,γT}(N)
        Nk = X0.Nk
        X1 = Array{G}(M)
        Xc = g0
        for k = 1:M
            X1[k] = Gconstr(args...)
        end
        Cc = Config(Nk, init=false)
        C1 = [Config(Nk, init=false) for k = 1:M]
        return new(N, Nk, X0, Xc, X1, Cc, C1)
    end
end

"""
    GraphLocalEntropy(N::Integer, M::Integer, γ::Float64, β::Float64, Gconstr, args...) <: DoubleGraph

A `DoubleGraph` which implements a "Local Entropy" model, given any other Ising model previously defined.
This simulates a replicated system in which each replica interacts with an extra "reference"
configuration.
This kind of graph can be simulated efficiently with [`rrrMC`](@ref).

`N` is the number of spins, `M` the number of replicas, `γ` the interaction strength,
`β` the inverse temperature. `GConstr` is the (original) graph constructor, and `args` the
arguments to the contructor.

This is similar to [`GraphRobustEnsemble`](@ref RRRMC.GraphRobustEnsemble), but the reference is simulated explicitly.

See also [`LEenergies`](@ref) and [`cenergy`](@ref).
"""
function GraphLocalEntropy(Nk::Integer, M::Integer, γ::Float64, β::Float64, Gconstr, args...)
    g0 = Gconstr(args...)
    G = typeof(g0)
    return GraphLocalEntropy{M,γ/β,G}(Nk * (M+1), g0, Gconstr, args...)
end

function update_cache!{M}(X::GraphLocalEntropy{M}, C::Config, move::Int)
    @extract X : X0 Xc X1 Cc C1
    k = mod1(move, M+1)
    i = (move - 1) ÷ (M+1) + 1

    if k == 1
        #spinflip!(Xc, Cc, i)
        spinflip!(Cc, i) # Xc cache not updated!
    else
        spinflip!(X1[k-1], C1[k-1], i)
    end

    update_cache!(X0, C, move)
end

function energy{M}(X::GraphLocalEntropy{M}, C::Config)
    # @assert X.N == C.N
    @extract X : Nk X0 Xc X1 Cc C1
    @extract C : s

    E = energy(X0, C)

    s1 = Cc.s
    for (i,j) = enumerate(1:(M+1):(1 + (M+1) * (Nk-1)))
        s1[i] = s[j]
    end
    energy(Xc, Cc)

    for k = 1:M
        s1 = C1[k].s
        for (i,j) = enumerate((k+1):(M+1):(k+1 + (M+1) * (Nk-1)))
            s1[i] = s[j]
        end
        E += energy(X1[k], C1[k])
    end

    return E
end

"""
    LEenergies(X::GraphLocalEntropy)

Returns a Vector with the individual energy (as defined by the original model)
of each replica in a [`GraphLocalEntropy`](@ref) graph.
"""
function LEenergies{M}(X::GraphLocalEntropy{M})
    @extract X : X1 C1

    Es = zeros(M)

    for k = 1:M
        Es[k] = energy(X1[k], C1[k])
    end

    return Es
end

function cenergy{M}(X::GraphLocalEntropy{M})
    @extract X : Xc Cc
    return energy(Xc, Cc)
end

function delta_energy_residual{M}(X::GraphLocalEntropy{M}, C::Config, move::Int)
    @extract X : X1 C1

    k = mod1(move, M+1)
    k == 1 && return 0.0

    i = (move - 1) ÷ (M+1) + 1

    return delta_energy(X1[k-1], C1[k-1], i)
end

function delta_energy(X::GraphLocalEntropy, C::Config, move::Int)
    return delta_energy(X.X0, C, move) +
           delta_energy_residual(X, C, move)
end

# This is type-unstable and inefficient. On the other hand, it is
# basically only written for testing purposes...
function neighbors{M}(X::GraphLocalEntropy{M}, i::Int)
    @extract X : X0 X1
    jts = neighbors(X0, i)

    k = mod1(i, M+1)

    k == 1 && return jts

    i1 = (i - 1) ÷ (M+1) + 1

    is = neighbors(X1[k-1], i1)

    return tuple(jts..., map(i->((i-1)*(M+1)+k), is)...)
end

function distances{M}(X::GraphLocalEntropy{M})
    @extract X : C1
    # @extract C : s
    # @assert Nk == C.N ÷ (M+1)

    # rs = BitVector[s[(k-1)*(M+1) + (2:(M+1):end)] for k = 1:M]

    # return [sum(s1 .⊻ s2) for s1 in rs, s2 in rs]
    return [sum(C1[k1].s .⊻ C1[k2].s) for k1 = 1:M, k2 = 1:M]
end

end # module
