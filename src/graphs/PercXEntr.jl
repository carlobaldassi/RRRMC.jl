# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module PercXEntr

using Random
using ExtractMacro
using ..Interface
using ..Common

using ..DeltaE.ArraySets

export GraphPercXEntr

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
    Δs::Vector{Int}
    newΔs::Vector{Int}
    ξsi::BitVector
    last_move::Int
    justΔ::Bool
    function Stabilities(P::Integer)
        Δs = zeros(P)
        newΔs = zeros(P)
        ξsi = BitArray(undef, P)
        return new(Δs, newΔs, ξsi, 0, false)
    end
end

function swap_Δs!(stab::Stabilities)
    stab.Δs, stab.newΔs = stab.newΔs, stab.Δs
    return stab
end

Δ2i(Δ, N) = (Δ + N) >>> 1 + 1

mutable struct GraphPercXEntr{λ} <: SimpleGraph{Float64}
    N::Int
    P::Int
    ξ::BitMatrix
    ξv::Vector{BitVector}
    stab::Stabilities
    Hs::Vector{Float64}
    tmps::BitVector
    function GraphPercXEntr{λ}(ξ::BitMatrix, ξv::Vector{BitVector}) where {λ}
        P, N = size(ξ)
        #TODO: check
        λ isa Float64 || throw(ArgumentError("λ must be a Float64, given a $(typeof(λ))"))
        isodd(N) || throw(ArgumentError("N must be odd, given: $N"))
        stab = Stabilities(P)
        Hs = [log(1 + exp(-2 * λ * Δ / √(N))) for Δ = -N:2:N]
        # for Δ = -N:2:N
        #     @assert Hs[Δ2i(Δ, N)] == log((1 + exp(-2 * λ * Δ / √(N))))
        # end
        tmps = BitVector(undef, N)
        return new{λ}(N, P, ξ, ξv, stab, Hs, tmps)
    end
end

"""
    GraphPercXEntr(N::Integer, P::Integer, λ::Float64) <: SimpleGraph{Float64}

A `SimpleGraph` implementing a single-layer binary perceptron with `N` binary (\$±1\$) synapses,
trained on `P` random i.i.d. \$±1\$ patterns.

The energy of the model is computed for each pattern using the cross-entropy loss (with a `λ`
parameter: ``\\log(1 + \\exp(-2λ σ⋅ξ / √N))``.

See also [`GraphPercStep`](@ref), [`GraphPercLinear`](@ref)
"""
GraphPercXEntr(N::Integer, P::Integer, λ::Float64) = GraphPercXEntr{λ}(gen_ξ(N, P)...)

GraphPercXEntr(X::GraphPercXEntr{λold}, λnew::Float64) where {λold} = GraphPercXEntr{λnew}(X.ξ, X.ξv)

function Base.empty!(stab::Stabilities)
    @extract stab : Δs
    fill!(Δs, 0)
    stab.last_move = 0
    stab.justΔ = false
    return stab
end

function energy(X::GraphPercXEntr{λ}, C::Config) where {λ}
    @assert X.N == C.N
    @extract C : s
    @extract X : N P ξv stab Hs tmps
    @extract stab : Δs

    E = 0.0
    empty!(stab)
    tmps = BitArray(undef, N)

    for a = 1:P
        # tmps[:] = ξ[a,:]
        # map!(⊻, tmps, s, tmps)
        # Δ = N - 2 * sum(tmps)
        map!(⊻, tmps, s, ξv[a])
        Δ = N - 2 * sum(tmps)
        Δs[a] = Δ
        # E += H(λ * Δ)
        E += Hs[Δ2i(Δ, N)]
    end

    return E
end

function update_cache!(X::GraphPercXEntr{λ}, C::Config, move::Int) where {λ}
    @assert X.N == C.N

    @extract X : stab
    @extract stab : last_move justΔ
    # NOTE: this only works if called right after delta_energy
    if justΔ && last_move == move
        swap_Δs!(stab)
        return
    end

    @extract C : s
    @extract X : N P ξ
    @extract stab : Δs ξsi
    @assert N == length(s)

    si = ~s[move]
    last_move ≠ move && unsafe_copyto!(ξsi, 1, ξ, (move-1)*P + 1, P)
    stab.last_move = move
    si && flipbits!(ξsi)
    # @assert ξsi == ξ[:,move] .⊻ si
    @inbounds for a = 1:P
        oldΔ = Δs[a]
        newΔ = oldΔ + (4 * ξsi[a] - 2)
        Δs[a] = newΔ
    end
    si && flipbits!(ξsi)
    stab.justΔ = false

    return
end

function delta_energy_naive(X::GraphPercXEntr, C::Config, move::Int)
    @assert C.N == X.N
    @assert 1 ≤ move ≤ N

    spinflip!(C, move)
    E1 = energy(X, C)
    spinflip!(C, move)
    E0 = energy(X, C)
    X.stab.justΔ = false
    return E1 - E0
end

function delta_energy(X::GraphPercXEntr{λ}, C::Config, move::Int) where {λ}
    @assert C.N == X.N
    @extract C : s
    @extract X : N P ξ stab Hs
    @extract stab : Δs newΔs ξsi last_move
    @assert 1 ≤ move ≤ N

    si = s[move]
    last_move ≠ move && unsafe_copyto!(ξsi, 1, ξ, (move-1)*P + 1, P)
    stab.last_move = move
    si && flipbits!(ξsi)
    # @assert ξsi == ξ[:,move] .⊻ si
    ΔE = 0.0
    @inbounds for a = 1:P
        oldΔ = Δs[a]
        newΔ = oldΔ + (4 * ξsi[a] - 2)
        newΔs[a] = newΔ
        # @assert -N ≤ newΔ ≤ N  (N, oldΔ, newΔ, ξsi[a])

        # ΔE += H(λ * newΔ) - H(λ * oldΔ)
        ΔE += Hs[Δ2i(newΔ, N)] - Hs[Δ2i(oldΔ, N)]
    end
    si && flipbits!(ξsi)
    # td = delta_energy_naive(X, s, stab, move)
    # ΔE == td || @show ΔE, td
    # @assert ΔE == td
    stab.justΔ = true
    return ΔE
end

function check_delta(X::GraphPercXEntr, C::Config, move::Int)
    delta = delta_energy(X, C, move)
    e0 = energy(X, C)
    spinflip!(X, C, move)
    e1 = energy(X, C)
    spinflip!(X, C, move)

    @assert isapprox(e1 - e0, delta, atol=1e-9) (e1, e0, delta, e1-e0)
end

function step_energy(X::GraphPercXEntr)
    @extract X : P stab
    @extract stab : Δs
    E = 0
    @inbounds for a = 1:P
        E += Δs[a] < 0
    end
    return E
end

neighbors(X::GraphPercXEntr, i::Int) = AllButOne(X.N, i)

end
