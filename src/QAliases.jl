# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module QAliases

using ..QT
using ..Empty
using ..SK
using ..EA
using ..Perc
using ..PercNaive

export GraphQ0T, GraphQSKT, GraphQSKNormalT, GraphQEAT, GraphQPercT, GraphQPercNaiveT

typealias GraphQ0T{fourK} GraphQuant{fourK,GraphEmpty}

"""
    GraphQ0T(N::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph

Shortcut for GraphQuant(N, M, Γ, β, GraphEmpty, N).

See [`GraphQuant`](@ref).

Intended for testing/debugging purposes.
"""
GraphQ0T(Nk::Integer, M::Integer, Γ::Float64, β::Float64) = GraphQuant(Nk, M, Γ, β, GraphEmpty, Nk)



typealias GraphQSKT{fourK} GraphQuant{fourK,GraphSK}

"""
    GraphQSKT(N::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph

Shortcut for GraphQuant(N, M, Γ, β, GraphSK, N), slightly optimized.

See [`GraphQuant`](@ref).
"""
GraphQSKT(Nk::Integer, M::Integer, Γ::Float64, β::Float64) = GraphQuant(Nk, M, Γ, β, GraphSK, SK.gen_J(Nk))


typealias GraphQSKNormalT{fourK} GraphQuant{fourK,GraphSKNormal}
GraphQSKNormalT(Nk::Integer, M::Integer, Γ::Float64, β::Float64) = GraphQuant(Nk, M, Γ, β, GraphSKNormal, SK.gen_J_gauss(Nk))



typealias GraphQEAT{fourK,twoD} GraphQuant{fourK,GraphEANormal{twoD}}

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

    GraphQuant(N, M, Γ, β, GraphEANormal{2D}, L, A, J)
end

function GraphQEAT(fname::AbstractString, M::Integer, Γ::Float64, β::Float64)
    L, D, A, J = EA.gen_AJ(fname)
    N = length(A)
    @assert N == L^D
    GraphQuant(N, M, Γ, β, GraphEANormal{2D}, L, A, J)
end

function GraphQEAT{twoD}(X::GraphEANormal{twoD}, M::Integer, Γ::Float64, β::Float64)
    @assert iseven(twoD)
    D = twoD ÷ 2
    N = X.N
    L = round(Int, N^(1/D))
    @assert L^D == N
    GraphQuant(N, M, Γ, β, GraphEANormal{twoD}, L, X.A, X.J)
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
