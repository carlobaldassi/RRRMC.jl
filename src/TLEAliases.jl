# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module TLEAliases

using Compat

using ..TLE
using ..Empty
using ..SK
using ..EA
using ..SAT

export Graph0TLE, GraphSKTLE, GraphEATLE, GraphSATTLE

@compat const Graph0TLE{M,γT,λT} = GraphTopologicalLocalEntropy{M,γT,λT,GraphEmpty}

# """
#     Graph0TLE(...)
#
# TODO
# """
Graph0TLE(Nk::Integer, M::Integer, γ::Float64, λ::Float64, β::Float64) = GraphTopologicalLocalEntropy(Nk, M, γ, λ, β, GraphEmpty, Nk)



@compat const GraphSKTLE{M,γT,λT} = GraphTopologicalLocalEntropy{M,γT,λT,GraphSK}

# """
#     GraphSKTLE(N::Integer, M::Integer, γ::Float64, λ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
GraphSKTLE(Nk::Integer, M::Integer, γ::Float64, λ::Float64, β::Float64) = GraphTopologicalLocalEntropy(Nk, M, γ, λ, β, GraphSK, SK.gen_J(Nk))




@compat const GraphEATLE{M,γT,λT,twoD} = GraphTopologicalLocalEntropy{M,γT,λT,GraphEANormal{twoD}}

# """
#     GraphEATLE(L::Integer, D::Integer, M::Integer, γ::Float64, λ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphEATLE(L::Integer, D::Integer, M::Integer, γ::Float64, λ::Float64, β::Float64)
    D ≥ 1 || throw(ArgumentError("D must be ≥ 0, given: $D"))
    A = EA.gen_EA(L, D)
    N = length(A)
    J = EA.gen_J(Float64, N, A) do
        4*rand() - 2
    end

    GraphTopologicalLocalEntropy(N, M, γ, λ, β, GraphEANormal{2D}, L, A, J)
end

function GraphEATLE(fname::AbstractString, M::Integer, γ::Float64, λ::Float64, β::Float64)
    L, D, A, J = EA.gen_AJ(fname)
    N = length(A)
    @assert N == L^D
    GraphTopologicalLocalEntropy(N, M, γ, λ, β, GraphEANormal{2D}, L, A, J)
end

function GraphEATLE{twoD}(X::GraphEANormal{twoD}, M::Integer, γ::Float64, λ::Float64, β::Float64)
    @assert iseven(twoD)
    D = twoD ÷ 2
    N = X.N
    L = round(Int, N^(1/D))
    @assert L^D == N
    GraphTopologicalLocalEntropy(N, M, γ, λ, β, GraphEANormal{twoD}, L, X.A, X.J)
end



@compat const GraphSATTLE{M,γT,λT} = GraphTopologicalLocalEntropy{M,γT,λT,GraphSAT}

# """
#     GraphSATTLE(N::Integer, K::Integer, α::Real, M::Integer, γ::Float64, λ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphSATTLE(N::Integer, K::Integer, α::Real, M::Integer, γ::Float64, λ::Float64, β::Float64)
    A, J = SAT.gen_randomKSAT(N, K, α)
    GraphTopologicalLocalEntropy(N, M, γ, λ, β, GraphSAT, N, A, J)
end

function GraphSATTLE(X::GraphSAT, M::Integer, γ::Float64, λ::Float64, β::Float64)
    N = X.N
    GraphTopologicalLocalEntropy(N, M, γ, λ, β, GraphSAT, N, X.A, X.J)
end

end # module
