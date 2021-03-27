# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module CommQu

using Random
using ExtractMacro
using DataStructures
using ..Interface
using ..Common

using ..DeltaE.ArraySets

export GraphCommQu

import ..Interface: energy, delta_energy, update_cache!, neighbors

function gen_ξ(N::Integer, P::Integer)
    ξ = BitArray(undef, P, N)
    ξv = Vector{BitVector}(undef, P)
    for a = 1:P
        ξv[a] = bitrand(N)
        ξ[a,:] = ξv[a]
    end
    # the ξv representation is for energy
    # the ξ representation is for delta_energy and update_cache!
    y = bitrand(P)
    return ξ, ξv, y
end

mutable struct Stabilities
    K2::Int
    p2::ArraySet
    m2::ArraySet
    Δ1s::Vector{Vector{Int}}
    Δ2s::Vector{Int}
    δΔ2Hsmin::Vector{MutableBinaryMinHeap{Float64}}
    δΔ2Hsmax::Vector{MutableBinaryMaxHeap{Float64}}
    ξsi::BitVector
    last_move::Int
    function Stabilities(K2::Integer, P::Integer)
        p2 = ArraySet(P)
        m2 = ArraySet(P)
        Δ1s = [zeros(P) for k in 1:K2]
        Δ2s = zeros(P)
        δΔ2Hsmin = [MutableBinaryMinHeap{Float64}() for a = 1:P]
        δΔ2Hsmax = [MutableBinaryMaxHeap{Float64}() for a = 1:P]
        ξsi = BitArray(undef, P)
        return new(K2, p2, m2, Δ1s, Δ2s, δΔ2Hsmin, δΔ2Hsmax, ξsi, 0)
    end
end

mutable struct GraphCommQu <: SimpleGraph{Int}
    N::Int
    K1::Int
    K2::Int
    P::Int
    ξ::BitMatrix
    ξv::Vector{BitVector}
    y::BitVector
    stab::Stabilities
    tmps::BitVector
    function GraphCommQu(K2::Int, ξ::BitMatrix, ξv::Vector{BitVector}, y::BitVector)
        P, N = size(ξ)
        @assert N % K2 == 0
        K1 = N ÷ K2
        @assert length(ξv) == P
        @assert all(length(ξ1) == N for ξ1 in ξv)
        @assert length(y) == P
        iseven(K1) || throw(ArgumentError("K1 must be even, given: $K1"))
        iseven(K2) || throw(ArgumentError("K2 must be even, given: $K2"))
        stab = Stabilities(K2, P)
        tmps = BitVector(undef, N)
        return new(N, K1, K2, P, ξ, ξv, y, stab, tmps)
    end
end

"""
    GraphCommQu(K1::Integer, K2::Integer, P::Integer; fc = false) <: SimpleGraph{Int}

A `SimpleGraph` implementing a two-layer binary committee machine with `K2` hidden units,
where each hidden unit has `K1` binary (\$±1\$) synapses, trained on `P` random i.i.d. \$±1\$
patterns. The output of each hidden unit is computed through a quadratic function \$x^2\$,
and then half of the `K2` units have a negative weight associated to their output.
Both `K1` and `K2` must be even.
The default is to generate a tree-like machine, but it can be made fully-connected with the
argument `fc=true`.

The energy of the model is computed as the number of misclassified patterns.
"""
function GraphCommQu(K1::Integer, K2::Integer, P::Integer; fc::Bool = false)
    Kin = fc ? K1 : K1 * K2
    ξ, ξv, y = gen_ξ(Kin, P)
    if fc
        ξ = repeat(ξ, outer=(1,K2))
        ξv = [repeat(ξ1, K2) for ξ1 in ξv]
    end
    GraphCommQu(K2, ξ, ξv, y)
end

function Base.empty!(stab::Stabilities)
    @extract stab : p2 m2 Δ1s Δ2s δΔ2Hsmin δΔ2Hsmax ξsi
    empty!(p2)
    empty!(m2)
    foreach(v->fill!(v, 0), Δ1s)
    fill!(Δ2s, 0)
    P = length(δΔ2Hsmin)
    for a = 1:P
        δΔ2Hsmin[a] = MutableBinaryMinHeap{Float64}()
        δΔ2Hsmax[a] = MutableBinaryMaxHeap{Float64}()
    end
    stab.last_move = 0
    return stab
end

function energy(X::GraphCommQu, C::Config)
    @assert X.N == C.N
    @extract C : s
    @extract X : N K1 K2 P ξv y stab tmps
    @extract stab : p2 m2 Δ1s Δ2s δΔ2Hsmin δΔ2Hsmax

    E = 0
    empty!(stab)
    tmps = BitArray(undef, K1)

    for a = 1:P
        o = 2y[a] - 1
        Δ2 = 0
        δΔ2Hmax = MutableBinaryMaxHeap{Float64}()
        δΔ2Hmin = MutableBinaryMinHeap{Float64}()
        for k = 1:K2
            c = o * (2 * (2k ≤ K2) - 1)
            unsafe_copyto!(tmps, 1, s, (k-1)*K1 + 1, K1)
            map!(⊻, tmps, tmps, ξv[a][(k-1)*K1 .+ (1:K1)])
            Δ1 = K1 - 2 * sum(tmps)
            Δ1s[k][a] = Δ1
            Δ2 += c * Δ1^2
            push!(δΔ2Hmax, 4 * (c + abs(c * Δ1)))
            push!(δΔ2Hmin, 4 * (c - abs(c * Δ1)))
        end
        Δ2s[a] = Δ2
        δΔ2Hsmax[a] = δΔ2Hmax
        δΔ2Hsmin[a] = δΔ2Hmin
        if Δ2 > 0
            Δ2 + first(δΔ2Hmin) ≤ 0 && push!(p2, a)
        else
            Δ2 + first(δΔ2Hmax) > 0 && push!(m2, a)
            E += 1
        end
    end

    return E
end

function update_cache!(X::GraphCommQu, C::Config, move::Int)
    @assert X.N == C.N
    @extract C : s
    @extract X : N K1 K2 P ξ y stab
    @extract stab : p2 m2 Δ1s Δ2s δΔ2Hsmin δΔ2Hsmax ξsi last_move
    @assert N == length(s)

    si = s[move]
    k2 = (move - 1) ÷ K1 + 1
    # k1 = (move - 1) % K1 + 1
    last_move ≠ move && unsafe_copyto!(ξsi, 1, ξ, (move-1)*P + 1, P)
    stab.last_move = move
    si && flipbits!(ξsi)
    Δ1k = Δ1s[k2]
    c = 2 * (2k2 ≤ K2) - 1
    # @assert ξsi == ξ[:,move] ⊻ si
    @inbounds for a = 1:P
        co = c * (2y[a] - 1)
        oldΔ2 = Δ2s[a]
        oldΔ1 = Δ1k[a]
        newΔ1 = oldΔ1 + (2 - 4 * ξsi[a])

        newΔ2 = oldΔ2 + co * (newΔ1^2 - oldΔ1^2)

        # oldδΔ2max = first(δΔ2Hsmax[a])
        # oldδΔ2min = first(δΔ2Hsmin[a])

        update!(δΔ2Hsmax[a], k2, 4 * (co + abs(co * newΔ1)))
        update!(δΔ2Hsmin[a], k2, 4 * (co - abs(co * newΔ1)))

        newδΔ2max = first(δΔ2Hsmax[a])
        newδΔ2min = first(δΔ2Hsmin[a])

        oldm2 = a ∈ m2
        newm2 = newΔ2 ≤ 0 && newΔ2 + newδΔ2max > 0
        oldp2 = a ∈ p2
        newp2 = newΔ2 > 0 && newΔ2 + newδΔ2min ≤ 0

        # oldm2 && (@assert oldΔ2 ≤ 0 && oldΔ2 + oldδΔ2max > 0)
        # oldp2 && (@assert oldΔ2 > 0 && oldΔ2 + oldδΔ2min ≤ 0)

        oldm2 && !newm2 && delete!(m2, a)
        !oldm2 && newm2 && push!(m2, a)
        oldp2 && !newp2 && delete!(p2, a)
        !oldp2 && newp2 && push!(p2, a)

        Δ1k[a] = newΔ1
        Δ2s[a] = newΔ2
    end
    si && flipbits!(ξsi)

    return
end

function delta_energy_naive(X::GraphCommQu, C::Config, move::Int)
    @assert C.N == X.N
    @extract C : s
    @extract X : N P ξ stab
    @assert 1 ≤ move ≤ N

    spinflip!(C, move)
    E1 = energy(X, C)
    spinflip!(C, move)
    E0 = energy(X, C)
    return E1 - E0
end

function delta_energy(X::GraphCommQu, C::Config, move::Int)
    @assert C.N == X.N
    @extract C : s
    @extract X : N K1 K2 P ξ y stab
    @extract stab : p2 m2 Δ1s Δ2s ξsi last_move
    @assert 1 ≤ move ≤ N

    si = s[move]
    k2 = (move - 1) ÷ K1 + 1
    # k1 = (move - 1) % K1 + 1
    c = 2 * (2k2 ≤ K2) - 1

    last_move ≠ move && unsafe_copyto!(ξsi, 1, ξ, (move-1)*P + 1, P)
    stab.last_move = move
    si && flipbits!(ξsi)
    Δ1k = Δ1s[k2]
    ΔE = 0
    @inbounds for a in p2
        co = c * (2y[a] - 1)
        if Δ2s[a] + co * 4 * (1 - (1 - 2*ξsi[a]) * Δ1k[a]) ≤ 0
            ΔE += 1
        end
    end
    @inbounds for a in m2
        co = c * (2y[a] - 1)
        if Δ2s[a] + co * 4 * (1 - (1 - 2*ξsi[a]) * Δ1k[a]) > 0
            ΔE -= 1
        end
    end
    si && flipbits!(ξsi)
    # td = delta_energy_naive(X, C, move)
    # ΔE == td || @show ΔE, td
    # @assert ΔE == td
    return ΔE
end

function check_delta(X::GraphCommQu, C::Config, move::Int)
    e0 = energy(X, C)
    delta = delta_energy(X, C, move)
    spinflip!(C, move)
    e1 = energy(X, C)
    spinflip!(C, move)

    (e1-e0) == delta || (@show e1,e0,delta,e1-e0; error())
end

neighbors(X::GraphCommQu, i::Int) = AllButOne(X.N, i)

end
