module QT

using ExtractMacro
using ..Interface
using ..Common

export GraphQT, Qenergy, transverse_mag

import ..Interface: energy, delta_energy, neighbors, allΔE

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
transverse_mag(X::DoubleGraph, C::Config, β::Float64) = transverse_mag(discr_graph(X), C, β)


const MAXDIGITS = 8

type GraphQT{fourK} <: DiscrGraph{Float64}
    N::Int
    M::Int
    Nk::Int
    function GraphQT(N::Integer, M::Integer)
        M > 2 || throw(ArgumentError("M must be greater than 2, given: $M"))
        isa(fourK, Float64) || throw(ArgumentError("invalid parameter fourK, expected Float64, given: $(typeof(fourK))"))
        round(fourK, MAXDIGITS) ≠ fourK && throw(ArgumentError("up to $MAXDIGITS decimal digits supported in fourK, given: $(fourK)"))
        N % M == 0 || throw(ArgumentError("N must be divisible by M, given: N=$N M=$M"))
        Nk = N ÷ M
        return new(N, M, Nk)
    end
end

@doc """
    GraphQT{fourK}(N::Integer, M::Integer) <: DiscrGraph

An auxiliary `DiscrGraph` used to implement the interactions in the
Suzuki-Trotter dimension when simulating quantum spin systems in a
transverse field.

It is only useful when implementing other graph types; see e.g. [`GraphQIsingT`](@ref).
""" -> GraphQT{fourK}(N::Integer, M::Integer)

GraphQT{oldK}(X::GraphQT{oldK}, newK::Float64) = GraphQT{newK}(X.N, X.M)

function energy0(X::GraphQT, C::Config)
    @assert X.N == C.N
    @extract X : M Nk
    @extract C : s
    n = 0
    for i = 1:Nk
        sj = s[i + (M-1) * Nk]
        for k = 1:M
            sk = s[i + (k-1) * Nk]
            n -= 1 - 2 * (sk $ sj)
            sj = sk
        end
    end
    return n
end

energy{fourK}(X::GraphQT{fourK}, C::Config) = energy0(X, C) * fourK / 4

function delta_energy{fourK}(X::GraphQT{fourK}, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : s

    k1, k2 = neighbors(X, move)
    sk = s[move]
    s1 = s[k1]
    s2 = s[k2]

    #Δ = ((2sk - 1) * ((2s1 - 1) + (2s2 - 1))) / 2
    #Δ = ((1 - 2 * (sk $ s1)) + (1 - 2 * (sk $ s2))) / 2
    Δ = ((sk $ ~s1) - (sk $ s2))
    #Δ = ((s1==sk) - (s2≠sk))
    #Δ = (s1 == s2) * (2 * (s1 == sk) - 1)

    return Δ * fourK
end

function neighbors(X::GraphQT, i::Int)
    @extract X : N Nk
    return (i - Nk + N * (i ≤ Nk), i + Nk - N * (i + Nk > N))
end
#neighbors(X::GraphQT, i::Int) = (mod1(i-X.Nk, X.N), mod1(i+X.Nk, X.N))

@generated allΔE{fourK}(::Type{GraphQT{fourK}}) = (0.0, fourK)

function transverse_mag{fourK}(X::GraphQT{fourK}, C::Config, β::Float64)
    @extract X : N M
    p = -energy0(X, C) / N

    # x = log(coth(β * Γ / M))
    x = β * fourK / 2

    return cosh(x) - p * sinh(x)
end

end
