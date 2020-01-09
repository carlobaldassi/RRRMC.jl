# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module CommStep

using Random
using ExtractMacro
using ..Interface
using ..Common

using ..DeltaE.ArraySets

export GraphCommStep

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
    return ξ, ξv
end

mutable struct Stabilities
    K2::Int
    p1::Vector{ArraySet}
    m1::Vector{ArraySet}
    p2::ArraySet
    m2::ArraySet
    Δ1s::Vector{Vector{Int}}
    Δ2s::Vector{Int}
    ξsi::BitVector
    last_move::Int
    function Stabilities(K2::Integer, P::Integer)
        p1 = [ArraySet(P) for k in 1:K2]
        m1 = [ArraySet(P) for k in 1:K2]
        p2 = ArraySet(P)
        m2 = ArraySet(P)
        Δ1s = [zeros(P) for k in 1:K2]
        Δ2s = zeros(P)
        ξsi = BitArray(undef, P)
        return new(K2, p1, m1, p2, m2, Δ1s, Δ2s, ξsi, 0)
    end
end

mutable struct GraphCommStep <: SimpleGraph{Int}
    N::Int
    K1::Int
    K2::Int
    P::Int
    ξ::BitMatrix
    ξv::Vector{BitVector}
    stab::Stabilities
    tmps::BitVector
    function GraphCommStep(K2::Int, ξ::BitMatrix, ξv::Vector{BitVector})
        P, N = size(ξ)
        @assert N % K2 == 0
        K1 = N ÷ K2
        @assert length(ξv) == P
        @assert all(length(ξ1) == N for ξ1 in ξv)
        isodd(K1) || throw(ArgumentError("K1 must be odd, given: $K1"))
        isodd(K2) || throw(ArgumentError("K2 must be odd, given: $K2"))
        stab = Stabilities(K2, P)
        tmps = BitVector(undef, N)
        return new(N, K1, K2, P, ξ, ξv, stab, tmps)
    end
end

"""
    GraphCommStep(K1::Integer, K2::Integer, P::Integer; fc = false) <: SimpleGraph{Int}

A `SimpleGraph` implementing a two-layer binary committee machine with `K2` hidden units,
where each hidden unit has `K1` binary (\$±1\$) synapses, trained on `P` random i.i.d. \$±1\$
patterns. The output of each hidden unit is computed through a step function \$sign(x)\$.
Both `K1` and `K2` must be odd.
The default is to generate a tree-like machine, but it can be made fully-connected with the
argument `fc=true`.

The energy of the model is computed as the number of misclassified patterns.
"""
function GraphCommStep(K1::Integer, K2::Integer, P::Integer; fc::Bool = false)
    Kin = fc ? K1 : K1 * K2
    ξ, ξv = gen_ξ(Kin, P)
    if fc
        ξ = repeat(ξ, outer=(1,K2))
        ξv = [repeat(ξ1, K2) for ξ1 in ξv]
    end
    GraphCommStep(K2, ξ, ξv)
end

function Base.empty!(stab::Stabilities)
    @extract stab : p1 m1 p2 m2 Δ1s Δ2s ξsi
    map!(empty!, p1, p1)
    map!(empty!, m1, m1)
    empty!(p2)
    empty!(m2)
    map!(v->fill!(v, 0), Δ1s, Δ1s)
    fill!(Δ2s, 0)
    stab.last_move = 0
    return stab
end

function energy(X::GraphCommStep, C::Config)
    @assert X.N == C.N
    @extract C : s
    @extract X : N K1 K2 P ξv stab tmps
    @extract stab : p1 m1 p2 m2 Δ1s Δ2s

    E = 0
    empty!(stab)
    tmps = BitArray(undef, K1)

    for a = 1:P
        Δ2 = 0
        for k = 1:K2
            unsafe_copyto!(tmps, 1, s, (k-1)*K1 + 1, K1)
            map!(⊻, tmps, tmps, ξv[a][(k-1)*K1 .+ (1:K1)])
            Δ1 = K1 - 2 * sum(tmps)
            Δ1s[k][a] = Δ1
            Δ2 += sign(Δ1)
            if Δ1 == 1
                push!(p1[k], a)
            elseif Δ1 == -1
                push!(m1[k], a)
            end
        end
        Δ2s[a] = Δ2
        if Δ2 == 1
            push!(p2, a)
        elseif Δ2 < 0
            Δ2 == -1 && push!(m2, a)
            E += 1
        end
    end

    return E
end

function update_cache!(X::GraphCommStep, C::Config, move::Int)
    @assert X.N == C.N
    @extract C : s
    @extract X : N K1 K2 P ξ stab
    @extract stab : p1 m1 p2 m2 Δ1s Δ2s ξsi last_move
    @assert N == length(s)

    si = s[move]
    k2 = (move - 1) ÷ K1 + 1
    # k1 = (move - 1) % K1 + 1
    last_move ≠ move && unsafe_copyto!(ξsi, 1, ξ, (move-1)*P + 1, P)
    stab.last_move = move
    si && flipbits!(ξsi)
    Δ1k = Δ1s[k2]
    p1k = p1[k2]
    m1k = m1[k2]
    # @assert ξsi == ξ[:,move] ⊻ si
    @inbounds for a = 1:P
        oldΔ2 = Δ2s[a]
        oldΔ1 = Δ1k[a]
        newΔ1 = oldΔ1 + (2 - 4 * ξsi[a])

        newΔ2 = oldΔ2
        if oldΔ1 > 1 && newΔ1 == 1
            push!(p1k, a)
        elseif oldΔ1 == 1
            delete!(p1k, a)
            newΔ1 == -1 && (push!(m1k, a); newΔ2 -= 2)
        elseif oldΔ1 == -1
            delete!(m1k, a)
            newΔ1 == 1 && (push!(p1k, a); newΔ2 += 2)
        elseif oldΔ1 < -1 && newΔ1 == -1
            push!(m1k, a)
        end
        Δ1k[a] = newΔ1

        newΔ2 == oldΔ2 && continue

        if oldΔ2 > 1 && newΔ2 == 1
            push!(p2, a)
        elseif oldΔ2 == 1
            delete!(p2, a)
            newΔ2 == -1 && push!(m2, a)
        elseif oldΔ2 == -1
            delete!(m2, a)
            newΔ2 == 1 && push!(p2, a)
        elseif oldΔ2 < -1 && newΔ2 == -1
            push!(m2, a)
        end
        Δ2s[a] = newΔ2
    end
    si && flipbits!(ξsi)

    return
end

function delta_energy_naive(X::GraphCommStep, C::Config, move::Int)
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

function delta_energy(X::GraphCommStep, C::Config, move::Int)
    @assert C.N == X.N
    @extract C : s
    @extract X : N K1 K2 P ξ stab
    @extract stab : p1 m1 p2 m2 ξsi last_move
    @assert 1 ≤ move ≤ N

    si = s[move]
    k2 = (move - 1) ÷ K1 + 1
    # k1 = (move - 1) % K1 + 1

    last_move ≠ move && unsafe_copyto!(ξsi, 1, ξ, (move-1)*P + 1, P)
    stab.last_move = move
    si && flipbits!(ξsi)
    ΔE = 0
    @inbounds for a in p2
        if a in p1[k2]
            ΔE += ~ξsi[a]
        end
    end
    @inbounds for a in m2
        if a in m1[k2]
            ΔE -= ξsi[a]
        end
    end
    si && flipbits!(ξsi)
    # td = delta_energy_naive(X, C, move)
    # ΔE == td || @show ΔE, td
    # @assert ΔE == td
    return ΔE
end

function check_delta(X::GraphCommStep, C::Config, move::Int)
    e0 = energy(X, C)
    delta = delta_energy(X, C, move)
    spinflip!(C, move)
    e1 = energy(X, C)
    spinflip!(C, move)

    (e1-e0) == delta || (@show e1,e0,delta,e1-e0; error())
end

neighbors(X::GraphCommStep, i::Int) = AllButOne(X.N, i)

end
