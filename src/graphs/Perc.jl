# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module Perc

using ExtractMacro
using Compat
using ..Interface
using ..Common

using ..DeltaE.ArraySets

export GraphPerc

import ..Interface: energy, delta_energy, update_cache!

function gen_ξ(N::Integer, P::Integer)
    ξ = BitArray(P, N)
    ξv = Vector{BitVector}(P)
    for a = 1:P
        ξv[a] = bitrand(N)
        ξ[a,:] = ξv[a]
    end
    # the ξv representation is for energy
    # the ξ representation is for delta_energy and update_cache!
    return ξ, ξv
end

type Stabilities
    p::ArraySet
    m::ArraySet
    Δs::Vector{Int}
    ξsi::BitVector
    last_move::Int
    function Stabilities(P::Integer)
        p = ArraySet(P)
        m = ArraySet(P)
        Δs = zeros(P)
        ξsi = BitArray(P)
        return new(p, m, Δs, ξsi, 0)
    end
end

type GraphPerc <: SimpleGraph{Int}
    N::Int
    P::Int
    ξ::BitMatrix
    ξv::Vector{BitVector}
    stab::Stabilities
    tmps::BitVector
    function GraphPerc(ξ::BitMatrix, ξv::Vector{BitVector})
        P, N = size(ξ)
        #TODO: check
        isodd(N) || throw(ArgumentError("N must be odd, given: $N"))
        stab = Stabilities(P)
        tmps = BitVector(N)
        return new(N, P, ξ, ξv, stab, tmps)
    end
end

"""
    GraphPerc(N::Integer, P::Integer) <: SimpleGraph{Int}

A `SimpleGraph` implementing a single-layer binary perceptron with `N` binary (\$±1\$) synapses,
trained on `P` random i.i.d. \$±1\$ patterns.

The energy of the model is computed for each pattern as the minimum number of weights which need
to be flipped in order to satisfy that pattern.

See also [`GraphPercNaive`](@ref)
"""
GraphPerc(N::Integer, P::Integer) = GraphPerc(gen_ξ(N, P)...)

function Base.empty!(stab::Stabilities)
    @extract stab : p m Δs ξsi
    empty!(p)
    empty!(m)
    fill!(Δs, 0)
    stab.last_move = 0
    return stab
end

function energy(X::GraphPerc, C::Config)
    @assert X.N == C.N
    @extract C : s
    @extract X : N P ξv stab tmps
    @extract stab : p m Δs

    E = 0
    empty!(stab)
    tmps = BitArray(N)

    for a = 1:P
        # tmps[:] = ξ[a,:]
        # map!(⊻, tmps, s, tmps)
        # Δ = N - 2 * sum(tmps)
        map!(⊻, tmps, s, ξv[a])
        Δ = N - 2 * sum(tmps)
        Δs[a] = Δ
        if Δ == 1
            push!(p, a)
        elseif Δ < 0
            push!(m, a)
            E += (-Δ-1) ÷ 2 + 1
        end
    end

    return E
end

function update_cache!(X::GraphPerc, C::Config, move::Int)
    @assert X.N == C.N
    @extract C : s
    @extract X : N P ξ stab
    @extract stab : p m Δs ξsi last_move
    @assert N == length(s)

    si = s[move]
    last_move ≠ move && unsafe_copy!(ξsi, 1, ξ, (move-1)*P + 1, P)
    stab.last_move = move
    si && flipbits!(ξsi)
    # @assert ξsi == ξ[:,move] ⊻ si
    @inbounds for a = 1:P
        oldΔ = Δs[a]
        newΔ = oldΔ + (2 - 4 * ξsi[a])

        if oldΔ > 1 && newΔ == 1
            push!(p, a)
        elseif oldΔ == 1
            delete!(p, a)
            newΔ < 0 && push!(m, a)
        elseif oldΔ < 0 && newΔ == 1
            delete!(m, a)
            push!(p, a)
        end
        Δs[a] = newΔ
    end
    si && flipbits!(ξsi)

    return
end

function delta_energy_naive(X::GraphPerc, C::Config, move::Int)
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

function delta_energy(X::GraphPerc, C::Config, move::Int)
    @assert C.N == X.N
    @extract C : s
    @extract X : N P ξ stab
    @extract stab : p m ξsi last_move
    @assert 1 ≤ move ≤ N

    si = s[move]
    last_move ≠ move && unsafe_copy!(ξsi, 1, ξ, (move-1)*P + 1, P)
    stab.last_move = move
    si && flipbits!(ξsi)
    # @assert ξsi == ξ[:,move] ⊻ si
    ΔE = 0
    @inbounds for a in p
        ΔE += ~ξsi[a]
    end
    @inbounds for a in m
        ΔE -= 2 * ξsi[a] - 1
    end
    si && flipbits!(ξsi)
    # td = delta_energy_naive(X, s, stab, move)
    # ΔE == td || @show ΔE, td
    # @assert ΔE == td
    return ΔE
end

function check_delta(X::GraphPerc, C::Config, move::Int)
    delta = delta_energy(X, C, move)
    e0 = energy(X, s)
    spinflip!(C, move)
    e1 = energy(X, s)
    spinflip!(C, move)

    (e1-e0) == delta || (@show e1,e0,delta,e1-e0; error())
end

end
