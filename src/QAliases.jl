# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module QAliases

using ..QT
using ..Empty
using ..SK
using ..EA
using ..PercLinear
using ..PercStep
using ..CommStep
using ..CommReLU
using ..CommQu

export GraphQ0T, GraphQSKT, GraphQSKNormalT, GraphQEAT,
       GraphQPercLinearT, GraphQPercStepT, GraphQCommStepT,
       GraphQCommReLUT, GraphQCommQuT

const GraphQ0T{fourK} = GraphQuant{fourK,GraphEmpty}

"""
    GraphQ0T(N::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph

Shortcut for GraphQuant(N, M, Γ, β, GraphEmpty, N).

See [`GraphQuant`](@ref).

Intended for testing/debugging purposes.
"""
GraphQ0T(Nk::Integer, M::Integer, Γ::Float64, β::Float64) = GraphQuant(Nk, M, Γ, β, GraphEmpty, Nk)



const GraphQSKT{fourK} = GraphQuant{fourK,GraphSK}

"""
    GraphQSKT(N::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph

Shortcut for GraphQuant(N, M, Γ, β, GraphSK, N), slightly optimized.

See [`GraphQuant`](@ref).
"""
GraphQSKT(Nk::Integer, M::Integer, Γ::Float64, β::Float64) = GraphQuant(Nk, M, Γ, β, GraphSK, SK.gen_J(Nk))


const GraphQSKNormalT{fourK} = GraphQuant{fourK,GraphSKNormal}
GraphQSKNormalT(Nk::Integer, M::Integer, Γ::Float64, β::Float64) = GraphQuant(Nk, M, Γ, β, GraphSKNormal, SK.gen_J_gauss(Nk))



const GraphQEAT{fourK,twoD} = GraphQuant{fourK,GraphEANormal{twoD}}

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

function GraphQEAT(X::GraphEANormal{twoD}, M::Integer, Γ::Float64, β::Float64) where {twoD}
    @assert iseven(twoD)
    D = twoD ÷ 2
    N = X.N
    L = round(Int, N^(1/D))
    @assert L^D == N
    GraphQuant(N, M, Γ, β, GraphEANormal{twoD}, L, X.A, X.J)
end

const GraphQPercLinearT{fourK} = GraphQuant{fourK,GraphPercLinear}

# """
#     GraphQPercLinearT(N::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphQPercLinearT(N::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64)
    ξ, ξv = PercLinear.gen_ξ(N, P)
    GraphQuant(N, M, Γ, β, GraphPercLinear, ξ, ξv)
end

function GraphQPercLinearT(X::GraphPercLinear, M::Integer, Γ::Float64, β::Float64)
    GraphQuant(X.N, M, Γ, β, GraphPercLinear, X.ξ, X.ξv)
end

const GraphQPercStepT{fourK} = GraphQuant{fourK,GraphPercStep}

# """
#     GraphQPercStepT(N::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphQPercStepT(N::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64)
    ξ, ξv = PercStep.gen_ξ(N, P)
    GraphQuant(N, M, Γ, β, GraphPercStep, ξ, ξv)
end

function GraphQPercStepT(X::GraphPercStep, M::Integer, Γ::Float64, β::Float64)
    GraphQuant(X.N, M, Γ, β, GraphPercStep, X.ξ, X.ξv)
end

const GraphQCommStepT{fourK} = GraphQuant{fourK,GraphCommStep}

# """
#     GraphQCommStepT(K1::Integer, K2::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphQCommStepT(K1::Integer, K2::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64; fc::Bool = false)
    N = K1 * K2
    Kin = fc ? K1 : N
    ξ, ξv = CommStep.gen_ξ(Kin, P)
    if fc
        ξ = repeat(ξ, outer=(1,K2))
        ξv = [repeat(ξ1, K2) for ξ1 in ξv]
    end
    GraphQuant(N, M, Γ, β, GraphCommStep, K2, ξ, ξv)
end

function GraphQCommStepT(X::GraphCommStep, M::Integer, Γ::Float64, β::Float64)
    GraphQuant(X.N, M, Γ, β, GraphCommStep, X.K2, X.ξ, X.ξv)
end

const GraphQCommReLUT{fourK} = GraphQuant{fourK,GraphCommReLU}

# """
#     GraphQCommReLUT(K1::Integer, K2::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphQCommReLUT(K1::Integer, K2::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64; fc::Bool = false)
    N = K1 * K2
    Kin = fc ? K1 : N
    ξ, ξv, y = CommReLU.gen_ξ(Kin, P)
    if fc
        ξ = repeat(ξ, outer=(1,K2))
        ξv = [repeat(ξ1, K2) for ξ1 in ξv]
    end
    GraphQuant(N, M, Γ, β, GraphCommReLU, K2, ξ, ξv, y)
end

function GraphQCommReLUT(X::GraphCommReLU, M::Integer, Γ::Float64, β::Float64)
    GraphQuant(X.N, M, Γ, β, GraphCommReLU, X.K2, X.ξ, X.ξv, X.y)
end

const GraphQCommQuT{fourK} = GraphQuant{fourK,GraphCommQu}

# """
#     GraphQCommQuT(K1::Integer, K2::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphQCommQuT(K1::Integer, K2::Integer, P::Integer, M::Integer, Γ::Float64, β::Float64; fc::Bool = false)
    N = K1 * K2
    Kin = fc ? K1 : N
    ξ, ξv, y = CommQu.gen_ξ(Kin, P)
    if fc
        ξ = repeat(ξ, outer=(1,K2))
        ξv = [repeat(ξ1, K2) for ξ1 in ξv]
    end
    GraphQuant(N, M, Γ, β, GraphCommQu, K2, ξ, ξv, y)
end

function GraphQCommQuT(X::GraphCommQu, M::Integer, Γ::Float64, β::Float64)
    GraphQuant(X.N, M, Γ, β, GraphCommQu, X.K2, X.ξ, X.ξv, X.y)
end

end # module
