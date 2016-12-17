# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module LEAliases

using ..LE
using ..Empty
using ..SK
using ..EA
using ..PercLinear
using ..PercStep

export Graph0LE, GraphSKLE, GraphEALE, GraphPercLinearLE, GraphPercStepLE

typealias Graph0LE{M,γT} GraphLocalEntropy{M,γT,GraphEmpty}

# """
#     Graph0LE(...)
#
# TODO
# """
Graph0LE(Nk::Integer, M::Integer, γ::Float64, β::Float64) = GraphLocalEntropy(Nk, M, γ, β, GraphEmpty, Nk)



typealias GraphSKLE{M,γT} GraphLocalEntropy{M,γT,GraphSK}

# """
#     GraphSKLE(N::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
GraphSKLE(Nk::Integer, M::Integer, γ::Float64, β::Float64) = GraphLocalEntropy(Nk, M, γ, β, GraphSK, SK.gen_J(Nk))




typealias GraphEALE{M,γT,twoD} GraphLocalEntropy{M,γT,GraphEANormal{twoD}}

# """
#     GraphEALE(L::Integer, D::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphEALE(L::Integer, D::Integer, M::Integer, γ::Float64, β::Float64)
    D ≥ 1 || throw(ArgumentError("D must be ≥ 0, given: $D"))
    A = EA.gen_EA(L, D)
    N = length(A)
    J = EA.gen_J(Float64, N, A) do
        4*rand() - 2
    end

    GraphLocalEntropy(N, M, γ, β, GraphEANormal{2D}, L, A, J)
end

function GraphEALE(fname::AbstractString, M::Integer, γ::Float64, β::Float64)
    L, D, A, J = EA.gen_AJ(fname)
    N = length(A)
    @assert N == L^D
    GraphLocalEntropy(N, M, γ, β, GraphEANormal{2D}, L, A, J)
end

function GraphEALE{twoD}(X::GraphEANormal{twoD}, M::Integer, γ::Float64, β::Float64)
    @assert iseven(twoD)
    D = twoD ÷ 2
    N = X.N
    L = round(Int, N^(1/D))
    @assert L^D == N
    GraphLocalEntropy(N, M, γ, β, GraphEANormal{twoD}, L, X.A, X.J)
end

typealias GraphPercLinearLE{M,γT} GraphLocalEntropy{M,γT,GraphPercLinear}

# """
#     GraphPercLinearLE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphPercLinearLE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64)
    ξ, ξv = PercLinear.gen_ξ(N, P)
    GraphLocalEntropy(N, M, γ, β, GraphPercLinear, ξ, ξv)
end

function GraphPercLinearLE(X::GraphPercLinear, M::Integer, γ::Float64, β::Float64)
    GraphLocalEntropy(X.N, M, γ, β, GraphPercLinear, X.ξ, X.ξv)
end

typealias GraphPercStepLE{M,γT} GraphLocalEntropy{M,γT,GraphPercStep}

# """
#     GraphPercStepLE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphPercStepLE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64)
    ξ, ξv = PercStep.gen_ξ(N, P)
    GraphLocalEntropy(N, M, γ, β, GraphPercStep, ξ, ξv)
end

function GraphPercStepLE(X::GraphPercStep, M::Integer, γ::Float64, β::Float64)
    GraphLocalEntropy(X.N, M, γ, β, GraphPercStep, X.ξ, X.ξv)
end

end # module
