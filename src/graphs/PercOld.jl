module PercOld

using ExtractMacro
using ..Interface
using ..Common
using ..QT

export GraphPercOld

import ..Interface: energy, delta_energy, update_cache!

function gen_ξ(N::Integer, P::Integer)
    ξ = BitVector[bitrand(N) for i = 1:P]
    return ξ
end

type Stabilities
    p::Vector{Int}
    m::Vector{Int}
    function Stabilities(N::Integer = 1)
        isodd(N) || throw(ArgumentError("N must be odd: $N"))
        s = round(Int, 5 * √N)
        p = Int[]
        m = Int[]
        sizehint!(p, s)
        sizehint!(m, N)
        return new(p, m)
    end
end

type GraphPercOld <: SimpleGraph{Int}
    N::Int
    P::Int
    ξ::Vector{BitVector}
    tmps::BitVector
    stab::Stabilities
    function GraphPercOld(ξ::Vector{BitVector}; check::Bool = true)
        P = length(ξ)
        @assert length(ξ) ≥ 1
        N = length(ξ[1])
        if check
            all(ξa->length(ξa) == N, ξ) || throw(ArgumentError("invalid ξ inner length, expected $N, given: $(unique(map(length,ξ)))"))
        end
        tmps = BitArray(N)
        stab = Stabilities(N)
        return new(N, P, ξ, tmps, stab)
    end
end

GraphPercOld(N::Integer, P::Integer) = GraphPercOld(gen_ξ(N, P), check=false)

function Base.empty!(stab::Stabilities)
    empty!(stab.p)
    empty!(stab.m)
    return stab
end

function energy(X::GraphPercOld, C::Config)
    @assert X.N == C.N
    @extract C : s
    @extract X : N P ξ tmps stab
    @extract stab : p m
    @assert N == length(s)

    E = 0
    empty!(stab)

    for a = 1:P
        ξa = ξ[a]
        map!($, tmps, s, ξa)
        Δ = N - 2 * sum(tmps)
        if Δ == 1
            push!(p, a)
        elseif Δ < 0
            push!(m, a)
            E += (-Δ-1) ÷ 2 + 1
        end
    end

    return E
end

# TODO: improve?
function update_cache!(X::GraphPercOld, C::Config, move::Int)
    energy(X, C)
end

function delta_energy_naive(X::GraphPercOld, C::Config, move::Int)
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

function delta_energy(X::GraphPercOld, C::Config, move::Int)
    @assert C.N == X.N
    @extract C : s
    @extract X : N P ξ stab
    @extract stab : p m
    @assert 1 ≤ move ≤ N

    si = s[move]
    Δ = 0
    @inbounds for a in p
        ξa = ξ[a]
        Δ += ξa[move] == si
    end
    @inbounds for a in m
        ξa = ξ[a]
        Δ -= 2 * (ξa[move] ≠ si) - 1
    end
    #td = delta_energy_naive(X, s, stab, move)
    #Δ == td || @show Δ, td
    #@assert Δ == td
    return Δ
end

function check_delta(X::GraphPercOld, C::Config, move::Int)
    delta = delta_energy(X, C, move)
    e0 = energy(X, s)
    spinflip!(C, move)
    e1 = energy(X, s)
    spinflip!(C, move)

    (e1-e0) == delta || (@show e1,e0,delta,e1-e0; error())
end

end
