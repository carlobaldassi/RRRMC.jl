module REAliases

using ..RE
using ..Empty
using ..IsingSK
using ..EA
using ..PercOld
using ..Perc
using ..PercNaive

export Graph0RE, GraphIsingRE, GraphEARE, GraphPercRE, GraphPercNaiveRE, GraphPercOldRE

typealias Graph0RE{M,γ,β} GraphRepl{M,γ,β,GraphEmpty}

# """
#     Graph0RE(...)
#
# TODO
# """
Graph0RE(Nk::Integer, M::Integer, γ::Float64, β::Float64) = GraphRepl(Nk, M, γ, β, GraphEmpty, Nk)



typealias GraphIsingRE{M,γ,β} GraphRepl{M,γ,β,GraphIsingSK}

# """
#     GraphIsingRE(N::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
GraphIsingRE(Nk::Integer, M::Integer, γ::Float64, β::Float64) = GraphRepl(Nk, M, γ, β, GraphIsingSK, IsingSK.gen_J(Nk))




typealias GraphEARE{M,γ,β,twoD} GraphRepl{M,γ,β,GraphEAContSimple{twoD}}

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

    GraphRepl(N, M, γ, β, GraphEAContSimple{2D}, L, A, J)
end

function GraphEARE(fname::AbstractString, M::Integer, γ::Float64, β::Float64)
    L, D, A, J = EA.gen_AJ(fname)
    N = length(A)
    @assert N == L^D
    GraphRepl(N, M, γ, β, GraphEAContSimple{2D}, L, A, J)
end

function GraphEARE{twoD}(X::GraphEAContSimple{twoD}, M::Integer, γ::Float64, β::Float64)
    @assert iseven(twoD)
    D = twoD ÷ 2
    N = X.N
    L = round(Int, N^(1/D))
    @assert L^D == N
    GraphRepl(N, M, γ, β, GraphEAContSimple{twoD}, L, X.A, X.J)
end

typealias GraphPercOldRE{M,γ,β} GraphRepl{M,γ,β,GraphPercOld}

# """
#     GraphPercOldRE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphPercOldRE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64)
    ξ = PercOld.gen_ξ(N, P)
    GraphRepl(N, M, γ, β, GraphPercOld, ξ)
end

function GraphPercOldRE(X::GraphPercOld, M::Integer, γ::Float64, β::Float64)
    GraphRepl(X.N, M, γ, β, GraphPercOld, X.ξ)
end

typealias GraphPercRE{M,γ,β} GraphRepl{M,γ,β,GraphPerc}

# """
#     GraphPercRE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphPercRE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64)
    ξ, ξv = Perc.gen_ξ(N, P)
    GraphRepl(N, M, γ, β, GraphPerc, ξ, ξv)
end

function GraphPercRE(X::GraphPerc, M::Integer, γ::Float64, β::Float64)
    GraphRepl(X.N, M, γ, β, GraphPerc, X.ξ, X.ξv)
end

typealias GraphPercNaiveRE{M,γ,β} GraphRepl{M,γ,β,GraphPercNaive}

# """
#     GraphPercNaiveRE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphPercNaiveRE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64)
    ξ, ξv = PercNaive.gen_ξ(N, P)
    GraphRepl(N, M, γ, β, GraphPercNaive, ξ, ξv)
end

function GraphPercNaiveRE(X::GraphPercNaive, M::Integer, γ::Float64, β::Float64)
    GraphRepl(X.N, M, γ, β, GraphPercNaive, X.ξ, X.ξv)
end

end # module
