# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module PercStep

using Random
using ExtractMacro
using ..Interface
using ..Common

using ..DeltaE.ArraySets

export GraphPercStep

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
    p::ArraySet
    m::ArraySet
    Δs::Vector{Int}
    ξsi::BitVector
    last_move::Int
    function Stabilities(P::Integer)
        p = ArraySet(P)
        m = ArraySet(P)
        Δs = zeros(P)
        ξsi = BitArray(undef, P)
        return new(p, m, Δs, ξsi, 0)
    end
end

mutable struct GraphPercStep <: SimpleGraph{Int}
    N::Int
    P::Int
    ξ::BitMatrix
    ξv::Vector{BitVector}
    stab::Stabilities
    tmps::BitVector
    function GraphPercStep(ξ::BitMatrix, ξv::Vector{BitVector})
        P, N = size(ξ)
        #TODO: check
        isodd(N) || throw(ArgumentError("N must be odd, given: $N"))
        stab = Stabilities(P)
        tmps = BitVector(undef, N)
        return new(N, P, ξ, ξv, stab, tmps)
    end
end

"""
    GraphPercStep(N::Integer, P::Integer) <: SimpleGraph{Int}

A `SimpleGraph` implementing a single-layer binary perceptron with `N` binary (\$±1\$) synapses,
trained on `P` random i.i.d. \$±1\$ patterns.

The energy of the model is computed as the number of misclassified patterns.

See also [`GraphPercLinear`](@ref).
"""
GraphPercStep(N::Integer, P::Integer) = GraphPercStep(gen_ξ(N, P)...)

function Base.empty!(stab::Stabilities)
    @extract stab : p m Δs ξsi
    empty!(p)
    empty!(m)
    fill!(Δs, 0)
    stab.last_move = 0
    return stab
end

function energy(X::GraphPercStep, C::Config)
    @assert X.N == C.N
    @extract C : s
    @extract X : N P ξv stab tmps
    @extract stab : p m Δs

    E = 0
    empty!(stab)
    tmps = BitArray(undef, N)

    for a = 1:P
        map!(⊻, tmps, s, ξv[a])
        Δ = N - 2 * sum(tmps)
        Δs[a] = Δ
        if Δ == 1
            push!(p, a)
        elseif Δ < 0
            Δ == -1 && push!(m, a)
            E += 1
        end
    end

    return E
end

function update_cache!(X::GraphPercStep, C::Config, move::Int)
    @assert X.N == C.N
    @extract C : s
    @extract X : N P ξ stab
    @extract stab : p m Δs ξsi last_move
    @assert N == length(s)

    si = s[move]
    last_move ≠ move && unsafe_copyto!(ξsi, 1, ξ, (move-1)*P + 1, P)
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
            newΔ == -1 && push!(m, a)
        elseif oldΔ == -1
            delete!(m, a)
            newΔ == 1 && push!(p, a)
        elseif oldΔ < -1 && newΔ == -1
            push!(m, a)
        end
        Δs[a] = newΔ
    end
    si && flipbits!(ξsi)

    return
end

function delta_energy_naive(X::GraphPercStep, C::Config, move::Int)
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

function delta_energy(X::GraphPercStep, C::Config, move::Int)
    @assert C.N == X.N
    @extract C : s
    @extract X : N P ξ stab
    @extract stab : p m ξsi last_move
    @assert 1 ≤ move ≤ N

    si = s[move]
    last_move ≠ move && unsafe_copyto!(ξsi, 1, ξ, (move-1)*P + 1, P)
    stab.last_move = move
    si && flipbits!(ξsi)
    # @assert ξsi == ξ[:,move] ⊻ si
    ΔE = 0
    @inbounds for a in p
        ΔE += ~ξsi[a]
    end
    @inbounds for a in m
        ΔE -= ξsi[a]
    end
    si && flipbits!(ξsi)
    # td = delta_energy_naive(X, s, stab, move)
    # ΔE == td || @show ΔE, td
    # @assert ΔE == td
    return ΔE
end

function check_delta(X::GraphPercStep, C::Config, move::Int)
    delta = delta_energy(X, C, move)
    e0 = energy(X, s)
    spinflip!(C, move)
    e1 = energy(X, s)
    spinflip!(C, move)

    (e1-e0) == delta || (@show e1,e0,delta,e1-e0; error())
end

neighbors(X::GraphPercStep, i::Int) = AllButOne(X.N, i)

end
