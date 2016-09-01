module QT

using ExtractMacro
using ..Interface
using ..Common

export GraphQT, Qenergy, transverse_mag,
       GraphQuant

import ..Interface: energy, delta_energy, neighbors, allΔE,
                    update_cache!, delta_energy_residual

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


# Add Transverse field to (almost) any AbstractGraph

type GraphQuant{fourK,G<:AbstractGraph} <: DoubleGraph{Float64}
    N::Int
    M::Int
    Nk::Int
    X0::GraphQT{fourK}
    X1::Vector{G}
    C1::Vector{Config}
    λ::Float64
    H0::Float64 # useless!!!
    β::Float64
    Γ::Float64
    function GraphQuant(N::Integer, M::Integer, H0::Real, β::Float64, Γ::Float64, g0::G, Gconstr, args...)
        X0 = GraphQT{fourK}(N, M)
        Nk = X0.Nk
        #J = gen_J(Nk)
        X1 = Array{G}(M)
        X1[1] = g0
        for k = 2:M
            X1[k] = Gconstr(args...)
        end
        C1 = [Config(Nk, init=false) for k = 1:M]
        λ = 1 / (M * √N)
        return new(N, M, Nk, X0, X1, C1, λ, H0, β, Γ)
    end
end

#  """
#      GraphQIsingT(N::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph
#
#  A `DoubleGraph` which implements a quantum Ising spin model in a transverse magnetic field,
#  using the Suzuki-Trotter transformation.
#  `N` is the number of spins, `M` the number of Suzuki-Trotter replicas, `Γ` the transverse
#  field, `β` the inverse temperature.
#  The graph is fully-connected, the interactions are random (\$J ∈ {-1,1}\$),
#  there are no external longitudinal fields.
#
#  See also [`Qenergy`](@ref).
#  """
function GraphQuant(Nk::Integer, M::Integer, Γ::Float64, β::Float64, Gconstr, args...)
    @assert Γ ≥ 0
    fourK = round(2/β * log(coth(β * Γ / M)), MAXDIGITS)
    H0 = Nk * M / 2β * log(sinh(2β * Γ / M) / 2) # useless!!!
    g0 = Gconstr(args...)
    G = typeof(g0)
    return GraphQuant{fourK,G}(Nk * M, M, H0, β, Γ, g0, Gconstr, args...)
end

function update_cache!(X::GraphQuant, C::Config, move::Int)
    @extract X : X1 Nk C1
    #@extract C : s
    k = (move - 1) ÷ Nk + 1
    i = mod1(move, Nk)

    spinflip!(X1[k], C1[k], i)

    #s1 = C1[k].s
    #copy!(s1, 1, s, (k-1)*Nk + 1, Nk)
    #@assert C1[k].s == s[((k-1)*Nk + 1):k*Nk]
end

function energy(X::GraphQuant, C::Config)
    @assert X.N == C.N
    @extract X : M Nk X0 X1 C1 λ H0
    @extract C : s

    E = energy(X0, C)

    for k = 1:M
        s1 = C1[k].s
        copy!(s1, 1, s, (k-1)*Nk + 1, Nk)
        E += energy(X1[k], C1[k]) * λ
    end

    return E - H0
end

function Qenergy(X::GraphQuant, C::Config)
    @assert X.N == C.N
    @extract X : M Nk X0 X1 C1 λ β Γ
    #@extract C : s

    E = -Γ * transverse_mag(X0, C, β)

    for k = 1:M
        #s1 = C1[k].s
        #copy!(s1, 1, s, (k-1)*Nk + 1, Nk)
        #@assert C1[k].s == s[((k-1)*Nk + 1):k*Nk]
        E += energy(X1[k], C1[k]) * λ / Nk
    end

    return E
end

function delta_energy_residual(X::GraphQuant, C::Config, move::Int)
    @extract X : Nk X1 C1 λ
    #@extract C : s

    k = (move - 1) ÷ Nk + 1
    #s1 = C1[k].s
    #copy!(s1, 1, s, (k-1)*Nk + 1, Nk)
    #@assert C1[k].s == s[((k-1)*Nk + 1):k*Nk]

    i = mod1(move, Nk)
    return delta_energy(X1[k], C1[k], i) * λ
end

function delta_energy(X::GraphQuant, C::Config, move::Int)
    return delta_energy(X.X0, C, move) +
           delta_energy_residual(X, C, move)
end

end
