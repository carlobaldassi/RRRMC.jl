# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module QT

using ExtractMacro
using ..Interface
using ..Common

using ...RRRMC # this is silly but it's required for correct cross-linking in docstrings, apparently

export GraphQT, Qenergy, transverse_mag,
       GraphQuant, Renergies, overlaps

import ..Interface: energy, delta_energy, neighbors, allΔE,
                    update_cache!, delta_energy_residual

import Base: iterate, length

"""
    Qenergy(X::DoubleGraph, C::Config)

When using the Suzuki-Trotter transformation to simulate quantum
systems in a transverse magnetic field with a replicated classical
system, this function should be used to obtain the average value of the
Hamiltonian observable (divided by the number of spins).
"""
Qenergy(::AbstractGraph, C::Config) = error("not implemented")

"""
    transverse_mag(X::DoubleGraph, C::Config, β::Float64)

When using the Suzuki-Trotter transformation to simulate quantum
systems in a transverse magnetic field with a replicated classical
system, this function should be used to obtain the average value of the
transverse magnetization observable.
"""
transverse_mag(X::DoubleGraph, C::Config, β::Float64) = transverse_mag(inner_graph(X), C, β)


const MAXDIGITS = 8

mutable struct GraphQT{fourK} <: DiscrGraph{Float64}
    N::Int
    M::Int
    Nk::Int
    function GraphQT{fourK}(N::Integer, M::Integer) where {fourK}
        M > 2 || throw(ArgumentError("M must be greater than 2, given: $M"))
        isa(fourK, Float64) || throw(ArgumentError("invalid parameter fourK, expected Float64, given: $(typeof(fourK))"))
        round(fourK, digits=MAXDIGITS) ≠ fourK && throw(ArgumentError("up to $MAXDIGITS decimal digits supported in fourK, given: $(fourK)"))
        N % M == 0 || throw(ArgumentError("N must be divisible by M, given: N=$N M=$M"))
        Nk = N ÷ M
        return new{fourK}(N, M, Nk)
    end
end

@doc """
    GraphQT{fourK}(N::Integer, M::Integer) <: DiscrGraph

An auxiliary `DiscrGraph` used to implement the interactions in the
Suzuki-Trotter dimension when simulating quantum spin systems in a
transverse field.

It is only useful when implementing other graph types; see [`GraphQuant`](@ref).
""" -> GraphQT{fourK}(N::Integer, M::Integer)

GraphQT(X::GraphQT{oldK}, newK::Float64) where {oldK} = GraphQT{newK}(X.N, X.M)

function energy0(X::GraphQT, C::Config)
    @assert X.N == C.N
    @extract X : M Nk
    @extract C : s
    n = 0
    for i = 1:Nk
        sj = s[i + (M-1) * Nk]
        for k = 1:M
            sk = s[i + (k-1) * Nk]
            n -= 1 - 2 * (sk ⊻ sj)
            sj = sk
        end
    end
    return n
end

energy(X::GraphQT{fourK}, C::Config) where {fourK} = energy0(X, C) * fourK / 4

function delta_energy(X::GraphQT{fourK}, C::Config, move::Int) where {fourK}
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : s

    k1, k2 = neighbors(X, move)
    sk = s[move]
    s1 = s[k1]
    s2 = s[k2]

    #Δ = ((2sk - 1) * ((2s1 - 1) + (2s2 - 1))) / 2
    #Δ = ((1 - 2 * (sk ⊻ s1)) + (1 - 2 * (sk ⊻ s2))) / 2
    Δ = ((sk ⊻ ~s1) - (sk ⊻ s2))
    #Δ = ((s1==sk) - (s2≠sk))
    #Δ = (s1 == s2) * (2 * (s1 == sk) - 1)

    return Δ * fourK
end

function neighbors(X::GraphQT, i::Int)
    @extract X : N Nk
    return (i - Nk + N * (i ≤ Nk), i + Nk - N * (i + Nk > N))
end
#neighbors(X::GraphQT, i::Int) = (mod1(i-X.Nk, X.N), mod1(i+X.Nk, X.N))

@generated allΔE(::Type{GraphQT{fourK}}) where {fourK} = (0.0, fourK)

function transverse_mag(X::GraphQT{fourK}, C::Config, β::Float64) where {fourK}
    @extract X : N M
    p = -energy0(X, C) / N

    # x = log(coth(β * Γ / M))
    x = β * fourK / 2

    return cosh(x) - p * sinh(x)
end


# Add Transverse field to (almost) any AbstractGraph

mutable struct GraphQuant{fourK,G<:AbstractGraph} <: DoubleGraph{DiscrGraph{Float64},Float64}
    N::Int
    M::Int
    Nk::Int
    X0::GraphQT{fourK}
    X1::Vector{G}
    C1::Vector{Config}
    β::Float64
    Γ::Float64
    function GraphQuant{fourK,G}(N::Integer, M::Integer, β::Float64, Γ::Float64, g0::G, Gconstr, args...) where {fourK,G}
        X0 = GraphQT{fourK}(N, M)
        Nk = X0.Nk
        #J = gen_J(Nk)
        X1 = Array{G}(undef, M)
        X1[1] = g0
        for k = 2:M
            X1[k] = Gconstr(args...)
        end
        C1 = [Config(Nk, init=false) for k = 1:M]
        return new{fourK,G}(N, M, Nk, X0, X1, C1, β, Γ)
    end
end

"""
    GraphQuant(N::Integer, M::Integer, Γ::Float64, β::Float64, Gconstr, args...) <: DoubleGraph

A `DoubleGraph` which implements a quantum Ising spin model in a transverse magnetic field,
using the Suzuki-Trotter transformation. This allows to model the quantum transverse field
case for any classical Ising model previously defined. This kind of graph can be simulated
efficiently with [`rrrMC`](@ref).

`N` is the number of spins, `M` the number of Suzuki-Trotter replicas, `Γ` the transverse
field, `β` the inverse temperature. `GConstr` is the (classical) graph constructor, and `args` the
arguments to the contructor.

See also [`Qenergy`](@ref) and [`transverse_mag`](@ref).
"""
function GraphQuant(Nk::Integer, M::Integer, Γ::Float64, β::Float64, Gconstr, args...)
    @assert Γ ≥ 0
    fourK = round(2/β * log(coth(β * Γ / M)), digits=MAXDIGITS)
    # H0 = Nk * M / 2β * log(sinh(2β * Γ / M) / 2) # useless!!!
    g0 = Gconstr(args...)
    G = typeof(g0)
    return GraphQuant{fourK,G}(Nk * M, M, β, Γ, g0, Gconstr, args...)
end

function update_cache!(X::GraphQuant, C::Config, move::Int)
    @extract X : X1 Nk C1
    #@extract C : s
    k = (move - 1) ÷ Nk + 1
    i = mod1(move, Nk)

    spinflip!(X1[k], C1[k], i)

    #s1 = C1[k].s
    #copyto!(s1, 1, s, (k-1)*Nk + 1, Nk)
    #@assert C1[k].s == s[((k-1)*Nk + 1):k*Nk]
end

function energy(X::GraphQuant, C::Config)
    @assert X.N == C.N
    @extract X : M Nk X0 X1 C1
    @extract C : s

    E = energy(X0, C)

    for k = 1:M
        s1 = C1[k].s
        copyto!(s1, 1, s, (k-1)*Nk + 1, Nk)
        E += energy(X1[k], C1[k]) / M
    end

    return E
end

function Renergies(X::GraphQuant)
    @extract X : M X1 C1

    Es = zeros(M)

    for k = 1:M
        Es[k] = energy(X1[k], C1[k])
    end

    return Es
end

function overlaps(X::GraphQuant)
    @extract X : M Nk C1

    aux = BitArray(undef, Nk)
    ovs = zeros(M ÷ 2)

    for k1 = 1:(M-1)
        @extract C1[k1] : s1=s
        for k2 = (k1+1):M
            @extract C1[k2] : s2=s
            map!(⊻, aux, s1, s2)
            δ = min(k2 - k1, M + k1 - k2)
            ovs[δ] += Nk - 2 * sum(aux)
        end
    end

    # for i = 1:Nk
    #     for k1 = 1:(M-1)
    #         s1 = s[i + (k1-1) * Nk]
    #         for k2 = (k1+1):M
    #             δ = min(k2 - k1, M + k1 - k2)
    #             s2 = s[i + (k2-1) * Nk]
    #             ovs[δ] += 1 - 2 * (sk1 ⊻ sk2)
    #         end
    #     end
    # end

    # note: for each node in the loop, there are
    # 2 edges linking to other nodes at a given distance δ,
    # except if M is even and δ=M÷2, in which case there is
    # only 1 edge. So the multiplicity is M for all overlaps
    # except for the last entry for even M, then it's M/2.
    for δ = 1:((M-1)÷2)
        ovs[δ] /= M * Nk
    end
    iseven(M) && (ovs[M÷2] /= M * Nk / 2)

    return ovs
end

function Qenergy(X::GraphQuant, C::Config)
    @assert X.N == C.N
    @extract X : M N X0 X1 C1 β Γ
    #@extract C : s

    E = -Γ * transverse_mag(X0, C, β)

    for k = 1:M
        #s1 = C1[k].s
        #copyto!(s1, 1, s, (k-1)*Nk + 1, Nk)
        #@assert C1[k].s == s[((k-1)*Nk + 1):k*Nk]
        E += energy(X1[k], C1[k]) / N
    end

    return E
end

function delta_energy_residual(X::GraphQuant, C::Config, move::Int)
    @extract X : M Nk X1 C1
    #@extract C : s

    k = (move - 1) ÷ Nk + 1
    #s1 = C1[k].s
    #copyto!(s1, 1, s, (k-1)*Nk + 1, Nk)
    #@assert C1[k].s == s[((k-1)*Nk + 1):k*Nk]

    i = mod1(move, Nk)
    return delta_energy(X1[k], C1[k], i) / M
end

function delta_energy(X::GraphQuant, C::Config, move::Int)
    return delta_energy(X.X0, C, move) +
           delta_energy_residual(X, C, move)
end

mutable struct QNeighbIter{T}
    tn1::Int
    tn2::Int
    off::Int
    r::T
    QNeighbIter(tn1, tn2, off, r::T) where {T} = new{T}(tn1, tn2, off, r)
end

iterate(qn::QNeighbIter) = (qn.tn1, (nothing, Val(1)))
iterate(qn::QNeighbIter, st::Tuple{Nothing,Val{1}}) = (qn.tn2, (nothing, Val(2)))
function iterate(qn::QNeighbIter, st::Tuple{Nothing,Val{2}})
    rr = iterate(qn.r)
    rr ≡ nothing && return nothing
    v, i = rr
    return (v + qn.off, (i, nothing))
end
function iterate(qn::QNeighbIter, st::Tuple{<:Any,Nothing})
    rr = iterate(qn.r, st[1])
    rr ≡ nothing && return nothing
    v, i = rr
    return (v + qn.off, (i, nothing))
end

length(qn::QNeighbIter) = 2 + length(qn.r)

function neighbors(X::GraphQuant, i::Int)
    @extract X : N Nk X0 X1
    j1, j2 = neighbors(X0, i)

    k = (i - 1) ÷ Nk + 1
    j = mod1(i, Nk)

    return QNeighbIter(j1, j2, (k-1)*Nk, neighbors(X1[k], j))
end

end # module
