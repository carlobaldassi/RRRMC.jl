module QAliases

using ..QT
using ..Empty
using ..IsingSK
using ..EA
using ..PercOld
using ..Perc
using ..PercNaive

export GraphQ0T, GraphQIsingT, GraphQIsingGaussT, GraphQEAT, GraphQPercT, GraphQPercNaiveT, GraphQPercOldT

typealias GraphQ0T{fourK} GraphQuant{fourK,GraphEmpty}

"""
    GraphQ0T(N::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph

A simple `DoubleGraph` which implements independent spins in a transverse magnetic field,
using the Suzuki-Trotter transformation.
`N` is the number of spins, `M` the number of Suzuki-Trotter replicas,
`Γ` the transverse field, `β` the inverse temperature.

Intended for testing/debugging purposes.

See also [`Qenergy`](@ref).
"""
GraphQ0T(Nk::Integer, M::Integer, Γ::Float64, β::Float64) = GraphQuant(Nk, M, Γ, β, GraphEmpty, Nk)



typealias GraphQIsingT{fourK} GraphQuant{fourK,GraphIsingSK}

"""
    GraphQIsingT(N::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph

A `DoubleGraph` which implements a quantum Ising spin model in a transverse magnetic field,
using the Suzuki-Trotter transformation.
`N` is the number of spins, `M` the number of Suzuki-Trotter replicas, `Γ` the transverse
field, `β` the inverse temperature.
The graph is fully-connected, the interactions are random (\$J ∈ {-1/√N,1/√N}\$),
there are no external longitudinal fields.

See also [`Qenergy`](@ref).
"""
GraphQIsingT(Nk::Integer, M::Integer, Γ::Float64, β::Float64) = GraphQuant(Nk, M, Γ, β, GraphIsingSK, IsingSK.gen_J(Nk))


typealias GraphQIsingGaussT{fourK} GraphQuant{fourK,GraphIsingSKGauss}
GraphQIsingGaussT(Nk::Integer, M::Integer, Γ::Float64, β::Float64) = GraphQuant(Nk, M, Γ, β, GraphIsingSKGauss, IsingSK.gen_J_gauss(Nk))



typealias GraphQEAT{fourK,twoD} GraphQuant{fourK,GraphEAContSimple{twoD}}

# """
#     GraphQEAT(L::Integer, D::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphQEAT(L::Integer, D::Integer, M::Integer, Γ::Float64, β::Float64)
    D ≥ 1 || throw(ArgumentError("D must be ≥ 0, given: $D"))
    A = EA.gen_EA(L, D)
    N = length(A)
    J = EA.gen_J(Float64, N, A) do
        4*rand() - 2
    end

    GraphQuant(N, M, Γ, β, GraphEAContSimple{2D}, L, A, J)
end

function GraphQEAT(fname::AbstractString, M::Integer, Γ::Float64, β::Float64)
    L, D, A, J = EA.gen_AJ(fname)
    N = length(A)
    @assert N == L^D
    GraphQuant(N, M, Γ, β, GraphEAContSimple{2D}, L, A, J)
end

function GraphQEAT{twoD}(X::GraphEAContSimple{twoD}, M::Integer, Γ::Float64, β::Float64)
    @assert iseven(twoD)
    D = twoD ÷ 2
    N = X.N
    L = round(Int, N^(1/D))
    @assert L^D == N
    GraphQuant(N, M, Γ, β, GraphEAContSimple{twoD}, L, X.A, X.J)
end

typealias GraphQPercOldT{fourK} GraphQuant{fourK,GraphPercOld}

# """
#     GraphQPercOldT(N::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphQPercOldT(N::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64)
    ξ = PercOld.gen_ξ(N, P)
    GraphQuant(N, M, Γ, β, GraphPercOld, ξ)
end

function GraphQPercOldT(X::GraphPercOld, M::Integer, Γ::Float64, β::Float64)
    GraphQuant(X.N, M, Γ, β, GraphPercOld, X.ξ)
end

typealias GraphQPercT{fourK} GraphQuant{fourK,GraphPerc}

# """
#     GraphQPercT(N::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphQPercT(N::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64)
    ξ, ξv = Perc.gen_ξ(N, P)
    GraphQuant(N, M, Γ, β, GraphPerc, ξ, ξv)
end

function GraphQPercT(X::GraphPerc, M::Integer, Γ::Float64, β::Float64)
    GraphQuant(X.N, M, Γ, β, GraphPerc, X.ξ, X.ξv)
end

typealias GraphQPercNaiveT{fourK} GraphQuant{fourK,GraphPercNaive}

# """
#     GraphQPercNaiveT(N::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphQPercNaiveT(N::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64)
    ξ, ξv = PercNaive.gen_ξ(N, P)
    GraphQuant(N, M, Γ, β, GraphPercNaive, ξ, ξv)
end

function GraphQPercNaiveT(X::GraphPercNaive, M::Integer, Γ::Float64, β::Float64)
    GraphQuant(X.N, M, Γ, β, GraphPercNaive, X.ξ, X.ξv)
end

end # module
