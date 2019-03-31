# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module TLE

using ExtractMacro
using ..Interface
using ..Common

using ...RRRMC # this is silly but it's required for correct cross-linking in docstrings, apparently

export GraphTLE, GraphTopologicalLocalEntropy, TLEenergies

import ..Interface: energy, delta_energy, neighbors, allΔE,
                    update_cache!, delta_energy_residual,
                    cenergy, distances

mutable struct GraphTLE{M,γT,λT} <: DiscrGraph{Float64}
    N::Int
    Nk::Int
    neighb::Vector{Vector{Int}}
    neighb_jcs::Vector{Vector{Int}}
    neighb_Us::Vector{Vector{UnitRange{Int}}}
    cache1::LocalFields{Int} # replica-center
    cache2::LocalFields{Int} # topological
    ncache::Vector{Int}

    function GraphTLE(N::Integer, neighb::Vector{Vector{Int}}) where {M,γT,λT}
        isa(M, Int) || throw(ArgumentError("invalid parameter M, expected Int, given: $(typeof(M))"))
        M > 2 || throw(ArgumentError("M must be greater than 2, given: $M"))
        isa(γT, Float64) || throw(ArgumentError("invalid parameter γT, expected Float64, given: $(typeof(γT))"))
        N % (M + 1) == 0 || throw(ArgumentError("N must be divisible by M+1, given: N=$N M=$M"))
        Nk = N ÷ (M+1)
        length(neighb) == Nk || throw(ArgumentError("invalid neighb argument, expected length $Nk, given: $(length(neighb))"))
        for (i,ni) in enumerate(neighb)
            all(j->(1 ≤ j ≤ Nk), ni) || throw(ArgumentError("invalid neighb[$i] argument, all elements must be between 1 and $Nk"))
            i ∈ ni && throw(ArgumentError("invalid neighb[$i], contains itself"))
            lni = length(ni)
            for i1 in 1:lni, i2 in (i1+1):lni
                ni[i1] ≠ ni[i2] || throw(ArgumentError("duplicated element $i1 in neighb[$i]"))
            end
        end
        for (i,ni) in enumerate(neighb)
            for i1 in ni
                i ∈ neighb[i1] || throw(ArgumentError("neighb is not symmetrical, $i1 ∈ ∂$i but $i ∉ ∂$i1"))
            end
        end

        neighb_jcs = [[1 + (i1-1) * (M+1) for i1 in ni] for ni in neighb]
        neighb_Us = [[jinterval(i1, M) for i1 in ni] for ni in neighb]

        cache1 = LocalFields{Int}(N)
        cache2 = LocalFields{Int}(N)
        ncache = empty!(Array{Int}(undef, (M+1) * maximum(length.(neighb))))
        return new{M,γT,λT}(N, Nk, neighb, neighb_jcs, neighb_Us, cache1, cache2, ncache)
    end
end

@doc """
    GraphTLE{M,γT,λT}(N::Integer) <: DiscrGraph

An auxiliary `DiscrGraph` used to implement the interactions in the
Local Entropy Ensemble.

It is only useful when implementing other graph types; see [`GraphTopologicalLocalEntropy`](@ref).
""" -> GraphTLE{M,γT,λT}(N::Integer)

GraphTLE(X::GraphTLE{M,oldγ,oldλ}, newγ::Float64, newλ::Float64) where {M,oldγ,oldλ} = GraphTLE{M,newγ,newλ}(X.N, X.neighb)

# @generated function ΔElist{M,γT}(::Type{GraphTLE{M,γT}})
#     Expr(:tuple, ntuple(d->fk(2*(d - 1 - (M-1) >>> 0x1) - iseven(M), γT), M)...)
# end
# lstind(μ̄::Int, M::Int) = (μ̄ + M-1) >>> 0x1 + 1
#
# function getk{M,γT}(X::GraphTLE{M,γT}, μ̄::Int)
#     @inbounds k = ΔElist(GraphTLE{M,γT})[lstind(μ̄, M)]
#     return k
# end

function energy(X::GraphTLE{M,γT,λT}, C::Config) where {M,γT,λT}
    # @assert X.N == C.N
    @extract X : Nk neighb cache1 cache2
    @extract cache1 : lfields1=lfields lfields_last1=lfields_last
    @extract cache2 : lfields2=lfields lfields_last2=lfields_last
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
            lfields1[j] = σc * σj
            μ += σj
        end
        f = σc * μ
        lfields1[jc] = f
        n -= f
    end
    cache1.move_last = 0
    fill!(lfields_last1, 0.0)

    r = 0
    fill!(lfields2, 0.0)
    for (i1, ni) in enumerate(neighb)
        jc1 = 1 + (i1-1) * (M+1)
        σc1 = 2s[jc1] - 1
        μ = 0
        for i2 in ni
            jc2 = 1 + (i2-1) * (M+1)
            σc2 = 2s[jc2] - 1
            σc = σc1 * σc2

            for k = 2:(M+1)
                j1 = k + (i1-1) * (M+1)
                j2 = k + (i2-1) * (M+1)
                σj = (2s[j1] - 1) * (2s[j2] - 1)
                lfields2[j1] += σc * σj
                lfields2[j2] += σc * σj
                μ += σc2 * σj
            end
        end
        f = σc1 * μ
        lfields2[jc1] = 2f
        r -= f
    end
    @assert r % 2 == 0
    r ÷= 2
    lfields2 .÷= 2
    # @show n,r
    # @show n * γT, r * λT
    cache2.move_last = 0
    fill!(lfields_last2, 0.0)

    return n * γT + r * λT
end

function kinterval(move::Integer, M::Integer)
    j0 = move - ((move-1) % (M+1))
    j1 = j0 + M
    return j0:j1
end

function jinterval(i::Integer, M::Integer)
    a = (M+1) * (i-1) + 1
    b = (M+1) * i
    return a:b
end

function update_cache!(X::GraphTLE{M,γT,λT}, C::Config, move::Int) where {M,γT,λT}
    # @assert X.N == C.N
    # @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : Nk neighb neighb_jcs neighb_Us cache1 cache2

    k = mod1(move, M+1)         # replica index
    i = (move-1) ÷ (M+1) + 1    # variable index
    jc = 1 + (i-1) * (M+1)      # corresponding reference
    if k ≠ 1
        jcs = neighb_jcs[i]
        # jcs = [1 + (i1-1) * (M+1) for i1 in neighb[i]]  # other neighboring references
        # js = [k + (i1-1) * (M+1) for i1 in neighb[i]]   # other neighboring variables
    else
        Ux = kinterval(move, M)                         # the whole slice
        # Us = [jinterval(i1, M) for i1 in neighb[i]]     # other neighboring slices
        Us = neighb_Us[i]
    end

    @extract cache1 : lfields1=lfields lfields_last1=lfields_last move_last1=move_last
    @extract cache2 : lfields2=lfields lfields_last2=lfields_last move_last2=move_last
    @assert move_last1 == move_last2
    if move_last1 == move
        @inbounds begin
            #for y in neighbors(X, move)
            if k ≠ 1
                lfields1[jc], lfields_last1[jc] = lfields_last1[jc], lfields1[jc]
                lfields1[move], lfields_last1[move] = lfields_last1[move], lfields1[move]
            else
                for y in Ux
                    lfields1[y], lfields_last1[y] = lfields_last1[y], lfields1[y]
                end
            end
            #lfields1[move] = -lfields1[move]
            #lfields_last1[move] = -lfields_last1[move]
        end

        @inbounds begin
            if k ≠ 1
                lfields2[jc], lfields_last2[jc] = lfields_last2[jc], lfields2[jc]
                lfields2[move], lfields_last2[move] = lfields_last2[move], lfields2[move]
                for jc1 in jcs
                    lfields2[jc1], lfields_last2[jc1] = lfields_last2[jc1], lfields2[jc1]
                    j1 = jc1 - 1 + k
                    lfields2[j1], lfields_last2[j1] = lfields_last2[j1], lfields2[j1]
                end
                # for j in js
                #     lfields2[j], lfields_last2[j] = lfields_last2[j], lfields2[j]
                # end
            else
                for y in Ux
                    lfields2[y], lfields_last2[y] = lfields_last2[y], lfields2[y]
                end
                for u in Us, y in u
                    lfields2[y], lfields_last2[y] = lfields_last2[y], lfields2[y]
                end
            end
        end
        return
    end

    @inbounds begin
        if k ≠ 1
            #jc = 1 + (i-1) * (M+1)
            σx = 2s[move] - 1
            σc = 2s[jc] - 1
            lfc = lfields1[jc]
            lfm = lfields1[move]
            lfields_last1[jc] = lfc
            lfields_last1[move] = lfm
            lfields1[jc] = lfc + 2 * (σc * σx)
            lfields1[move] = -lfm
        else
            for y in Ux
                lfy = lfields1[y]
                lfields_last1[y] = lfy
                lfields1[y] = -lfy
            end
        end
        #lfm = lfields1[move]
        #lfields_last1[move] = lfm
        #lfields1[move] = -lfm
    end
    cache1.move_last = move

    # lfields_bk1 = copy(lfields1)
    # energy(X, C)
    # lfields_bk1 ≠ lfields1 && @show move hcat(lfields1,lfields_bk1) find(lfields1 .≠ lfields_bk1)
    # @assert lfields_bk1 == lfields1

    @inbounds begin
        if k ≠ 1
            σx = 2s[move] - 1
            lfm = lfields2[move]
            lfields_last2[move] = lfm
            lfields2[move] = -lfm

            σc = 2s[jc] - 1
            lfc = lfields2[jc]
            lfields_last2[jc] = lfc
            # for (jc1,j1) in zip(jcs, js)
            for jc1 in jcs
                j1 = jc1 - 1 + k
                lfj1 = lfields2[j1]
                lfc1 = lfields2[jc1]
                lfields_last2[j1] = lfj1
                lfields_last2[jc1] = lfc1
                σj1 = 2s[j1] - 1
                σc1 = 2s[jc1] - 1
                σ = σx * σc * σj1 * σc1
                lfields2[j1] = lfj1 + 2 * σ
                lfields2[jc1] = lfc1 + 2 * σ
                lfc += 2σ
            end
            lfields2[jc] = lfc
        else
            for y in Ux
                lfy = lfields2[y]
                lfields_last2[y] = lfy
                lfields2[y] = -lfy
            end
            σc = 2s[jc] - 1
            for Ux1 in Us
                jc1 = first(Ux1)
                lfc1 = lfields2[jc1]
                lfields_last2[jc1] = lfc1
                σc1 = 2s[jc1] - 1
                for j1 = (jc1+1):last(Ux1)
                    lfj1 = lfields2[j1]
                    lfields_last2[j1] = lfj1
                    σj1 = 2s[j1] - 1
                    σx = 2s[jc+(j1-jc1)] - 1
                    σ = σx * σc * σj1 * σc1
                    lfields2[j1] = lfj1 + 2 * σ
                    lfc1 += 2 * σ
                end
                lfields2[jc1] = lfc1
            end
        end
    end
    cache2.move_last = move

    return
end

@inline function delta_energy(X::GraphTLE{M,γT,λT}, C::Config, move::Int) where {M,γT,λT}
    # @assert X.N == C.N
    # @assert 1 ≤ move ≤ C.N
    @extract X : cache1 cache2
    @extract cache1 : lfields1=lfields
    @extract cache2 : lfields2=lfields

    # @show 2γT * lfields1[move], 2λT * lfields2[move]
    @inbounds Δ = 2γT * lfields1[move] + 2λT * lfields2[move]
    return Δ
end

@inline function neighbors(X::GraphTLE{M}, j::Int) where {M}
    @extract X : neighb ncache
    k = mod1(j, M+1)         # replica index
    i = (j-1) ÷ (M+1) + 1    # variable index
    if k == 1
        resize!(ncache, M + (M+1) * length(neighb[i]))
        ncache[1:M] = j + (1:M)
        for (n,i1) in enumerate(neighb[i])
            ncache[M + jinterval(n, M)] = jinterval(i1, M)
        end
        # return (j+1):(j+M)
    else
        resize!(ncache, 1 + 2 * length(neighb[i]))
        ncache[1] = j - k + 1
        for (n,i1) in enumerate(neighb[i])
            ncache[2n] = 1 + (i1-1) * (M+1)
            ncache[2n+1] = k + (i1-1) * (M+1)
        end
    end
    return ncache
end

function allΔE(X::GraphTLE{M,γT,λT}) where {M,γT,λT}
    @extract X : neighb
    ΔE1 = iseven(M) ? insert!([4*d*γT for d = 0:(M÷2)], 2, 2γT) :
                      [2*(2d-1)*γT for d = 1:((M+1)÷2)]

    mn = M * maximum(length.(neighb))
    ΔE2 = [2*d*λT for d = -mn:mn]
    # @show ΔE1
    # @show ΔE2

    ΔE = sort!(unique(abs.(ΔE1 .+ ΔE2')))
    return tuple(ΔE...)
end

# Replicate an existsing graph

mutable struct GraphTopologicalLocalEntropy{M,γT,λT,G<:AbstractGraph} <: DoubleGraph{DiscrGraph{Float64},Float64}
    N::Int
    Nk::Int
    X0::GraphTLE{M,γT,λT}
    Xc::G
    X1::Vector{G}
    Cc::Config
    C1::Vector{Config}
    function GraphTopologicalLocalEntropy(N::Integer, neighb::Vector{Vector{Int}}, g0::G, Gconstr, args...) where {M,γT,λT,G}
        X0 = GraphTLE{M,γT,λT}(N, neighb)
        Nk = X0.Nk
        X1 = Array{G}(undef, M)
        Xc = g0
        for k = 1:M
            X1[k] = Gconstr(args...)
        end
        Cc = Config(Nk, init=false)
        C1 = [Config(Nk, init=false) for k = 1:M]
        return new{M,γT,λT,G}(N, Nk, X0, Xc, X1, Cc, C1)
    end
end

"""
    GraphTopologicalLocalEntropy(N::Integer, M::Integer, γ::Float64, λ::Float64, β::Float64, Gconstr, args...) <: DoubleGraph

A `DoubleGraph` which implements a "Local Entropy" model, given any other Ising model previously defined.
This simulates a replicated system in which each replica interacts with an extra "reference"
configuration.
This kind of graph can be simulated efficiently with [`rrrMC`](@ref).

`N` is the number of spins, `M` the number of replicas, `γ` the interaction strength of the replicas with the reference,
`λ` is the topological interaction strenght,
`β` the inverse temperature. `GConstr` is the (original) graph constructor, and `args` the
arguments to the contructor.

This is similar to [`GraphLocalEntropy`](@ref RRRMC.GraphLocalEntropy), but with the additional topological interaction.

See also [`TLEenergies`](@ref) and [`cenergy`](@ref).
"""
function GraphTopologicalLocalEntropy(Nk::Integer, M::Integer, γ::Float64, λ::Float64, β::Float64, Gconstr, args...)
    g0 = Gconstr(args...)
    getN(g0) == Nk || throw(ArgumentError("expected $Nk spins, got $(getN(g0))"))
    G = typeof(g0)
    neighb = [Int[neighbors(g0, i)...] for i = 1:Nk]
    return GraphTopologicalLocalEntropy{M,γ/β,λ/β,G}(Nk * (M+1), neighb, g0, Gconstr, args...)
end

function update_cache!(X::GraphTopologicalLocalEntropy{M}, C::Config, move::Int) where {M}
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

function energy(X::GraphTopologicalLocalEntropy{M}, C::Config) where {M}
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
    TLEenergies(X::GraphTopologicalLocalEntropy)

Returns a Vector with the individual energy (as defined by the original model)
of each replica in a [`GraphTopologicalLocalEntropy`](@ref) graph.
"""
function TLEenergies(X::GraphTopologicalLocalEntropy{M}) where {M}
    @extract X : X1 C1

    Es = zeros(M)

    for k = 1:M
        Es[k] = energy(X1[k], C1[k])
    end

    return Es
end

function cenergy(X::GraphTopologicalLocalEntropy{M}) where {M}
    @extract X : Xc Cc
    return energy(Xc, Cc)
end

function delta_energy_residual(X::GraphTopologicalLocalEntropy{M}, C::Config, move::Int) where {M}
    @extract X : X1 C1

    k = mod1(move, M+1)
    k == 1 && return 0.0

    i = (move - 1) ÷ (M+1) + 1

    return delta_energy(X1[k-1], C1[k-1], i)
end

function delta_energy(X::GraphTopologicalLocalEntropy, C::Config, move::Int)
    return delta_energy(X.X0, C, move) +
           delta_energy_residual(X, C, move)
end

# This is type-unstable and inefficient. On the other hand, it is
# basically only written for testing purposes...
function neighbors(X::GraphTopologicalLocalEntropy{M}, i::Int) where {M}
    @extract X : X0 X1
    jts = neighbors(X0, i)

    k = mod1(i, M+1)

    k == 1 && return jts

    i1 = (i - 1) ÷ (M+1) + 1

    is = neighbors(X1[k-1], i1)

    return tuple(jts..., map(i->((i-1)*(M+1)+k), is)...)
end

function distances(X::GraphTopologicalLocalEntropy{M}) where {M}
    @extract X : C1
    # @extract C : s
    # @assert Nk == C.N ÷ (M+1)

    # rs = BitVector[s[(k-1)*(M+1) + (2:(M+1):end)] for k = 1:M]

    # return [sum(s1 .⊻ s2) for s1 in rs, s2 in rs]
    return [sum(C1[k1].s .⊻ C1[k2].s) for k1 = 1:M, k2 = 1:M]
end

end # module
