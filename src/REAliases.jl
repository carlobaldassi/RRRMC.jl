# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module REAliases

using ..RE
using ..Empty
using ..SK
using ..EA
using ..SAT
using ..PercLinear
using ..PercStep
using ..CommStep

export Graph0RE, GraphSKRE, GraphEARE, GraphSATRE, GraphPercLinearRE, GraphPercStepRE, GraphCommStepRE

const Graph0RE{M,γ,β} = GraphRobustEnsemble{M,γ,β,GraphEmpty}

# """
#     Graph0RE(...)
#
# TODO
# """
Graph0RE(Nk::Integer, M::Integer, γ::Float64, β::Float64) = GraphRobustEnsemble(Nk, M, γ, β, GraphEmpty, Nk)



const GraphSKRE{M,γ,β} = GraphRobustEnsemble{M,γ,β,GraphSK}

# """
#     GraphSKRE(N::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
GraphSKRE(Nk::Integer, M::Integer, γ::Float64, β::Float64) = GraphRobustEnsemble(Nk, M, γ, β, GraphSK, SK.gen_J(Nk))




const GraphEARE{M,γ,β,twoD} = GraphRobustEnsemble{M,γ,β,GraphEANormal{twoD}}

# """
#     GraphEARE(L::Integer, D::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphEARE(L::Integer, D::Integer, M::Integer, γ::Float64, β::Float64)
    D ≥ 1 || throw(ArgumentError("D must be ≥ 0, given: $D"))
    A = EA.gen_EA(L, D)
    N = length(A)
    J = EA.gen_J(Float64, N, A) do
        4*rand() - 2
    end

    GraphRobustEnsemble(N, M, γ, β, GraphEANormal{2D}, L, A, J)
end

function GraphEARE(fname::AbstractString, M::Integer, γ::Float64, β::Float64)
    L, D, A, J = EA.gen_AJ(fname)
    N = length(A)
    @assert N == L^D
    GraphRobustEnsemble(N, M, γ, β, GraphEANormal{2D}, L, A, J)
end

function GraphEARE(X::GraphEANormal{twoD}, M::Integer, γ::Float64, β::Float64) where {twoD}
    @assert iseven(twoD)
    D = twoD ÷ 2
    N = X.N
    L = round(Int, N^(1/D))
    @assert L^D == N
    GraphRobustEnsemble(N, M, γ, β, GraphEANormal{twoD}, L, X.A, X.J)
end

const GraphSATRE{M,γ,β} = GraphRobustEnsemble{M,γ,β,GraphSAT}

# """
#     GraphSATRE(N::Integer, K::Integer, α::Real, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphSATRE(N::Integer, K::Integer, α::Real, M::Integer, γ::Float64, β::Float64)
    A, J = SAT.gen_randomKSAT(N, K, α)
    GraphRobustEnsemble(N, M, γ, β, GraphSAT, N, A, J)
end

function GraphSATRE(X::GraphSAT, M::Integer, γ::Float64, β::Float64)
    N = X.N
    GraphRobustEnsemble(N, M, γ, β, GraphSAT, N, X.A, X.J)
end

const GraphPercLinearRE{M,γ,β} = GraphRobustEnsemble{M,γ,β,GraphPercLinear}

# """
#     GraphPercLinearRE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphPercLinearRE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64)
    ξ, ξv = PercLinear.gen_ξ(N, P)
    GraphRobustEnsemble(N, M, γ, β, GraphPercLinear, ξ, ξv)
end

function GraphPercLinearRE(X::GraphPercLinear, M::Integer, γ::Float64, β::Float64)
    GraphRobustEnsemble(X.N, M, γ, β, GraphPercLinear, X.ξ, X.ξv)
end

const GraphPercStepRE{M,γ,β} = GraphRobustEnsemble{M,γ,β,GraphPercStep}

# """
#     GraphPercStepRE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphPercStepRE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64)
    ξ, ξv = PercStep.gen_ξ(N, P)
    GraphRobustEnsemble(N, M, γ, β, GraphPercStep, ξ, ξv)
end

function GraphPercStepRE(X::GraphPercStep, M::Integer, γ::Float64, β::Float64)
    GraphRobustEnsemble(X.N, M, γ, β, GraphPercStep, X.ξ, X.ξv)
end

const GraphCommStepRE{M,γ,β} = GraphRobustEnsemble{M,γ,β,GraphCommStep}

# """
#     GraphCommStepRE(K1::Integer, K2::Integer, P::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphCommStepRE(K1::Integer, K2::Integer, P::Integer, M::Integer, γ::Float64, β::Float64)
    ξ, ξv = CommStep.gen_ξ(K1, P)
    N = K1 * K2
    GraphRobustEnsemble(N, M, γ, β, GraphCommStep, K2, ξ, ξv)
end

function GraphCommStepRE(X::GraphCommStep, M::Integer, γ::Float64, β::Float64)
    GraphRobustEnsemble(X.N, M, γ, β, GraphCommStep, X.K2, X.ξ, X.ξv)
end

end # module
