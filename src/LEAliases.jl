module LEAliases

using ..LE
using ..Empty
using ..IsingSK
using ..EA
using ..PercOld
using ..Perc
using ..PercNaive

export Graph0LE, GraphIsingLE, GraphEALE, GraphPercLE, GraphPercNaiveLE, GraphPercOldLE

typealias Graph0LE{M,γT} GraphLocEntr{M,γT,GraphEmpty}

# """
#     Graph0LE(...)
#
# TODO
# """
Graph0LE(Nk::Integer, M::Integer, γ::Float64, β::Float64) = GraphLocEntr(Nk, M, γ, β, GraphEmpty, Nk)



typealias GraphIsingLE{M,γT} GraphLocEntr{M,γT,GraphIsingSK}

# """
#     GraphIsingLE(N::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
GraphIsingLE(Nk::Integer, M::Integer, γ::Float64, β::Float64) = GraphLocEntr(Nk, M, γ, β, GraphIsingSK, IsingSK.gen_J(Nk))




typealias GraphEALE{M,γT,twoD} GraphLocEntr{M,γT,GraphEAContSimple{twoD}}

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

    GraphLocEntr(N, M, γ, β, GraphEAContSimple{2D}, L, A, J)
end

function GraphEALE(fname::AbstractString, M::Integer, γ::Float64, β::Float64)
    L, D, A, J = EA.gen_AJ(fname)
    N = length(A)
    @assert N == L^D
    GraphLocEntr(N, M, γ, β, GraphEAContSimple{2D}, L, A, J)
end

function GraphEALE{twoD}(X::GraphEAContSimple{twoD}, M::Integer, γ::Float64, β::Float64)
    @assert iseven(twoD)
    D = twoD ÷ 2
    N = X.N
    L = round(Int, N^(1/D))
    @assert L^D == N
    GraphLocEntr(N, M, γ, β, GraphEAContSimple{twoD}, L, X.A, X.J)
end

typealias GraphPercOldLE{M,γT} GraphLocEntr{M,γT,GraphPercOld}

# """
#     GraphPercOldLE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphPercOldLE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64)
    ξ = PercOld.gen_ξ(N, P)
    GraphLocEntr(N, M, γ, β, GraphPercOld, ξ)
end

function GraphPercOldLE(X::GraphPercOld, M::Integer, γ::Float64, β::Float64)
    GraphLocEntr(X.N, M, γ, β, GraphPercOld, X.ξ)
end

typealias GraphPercLE{M,γT} GraphLocEntr{M,γT,GraphPerc}

# """
#     GraphPercLE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphPercLE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64)
    ξ, ξv = Perc.gen_ξ(N, P)
    GraphLocEntr(N, M, γ, β, GraphPerc, ξ, ξv)
end

function GraphPercLE(X::GraphPerc, M::Integer, γ::Float64, β::Float64)
    GraphLocEntr(X.N, M, γ, β, GraphPerc, X.ξ, X.ξv)
end

typealias GraphPercNaiveLE{M,γT} GraphLocEntr{M,γT,GraphPercNaive}

# """
#     GraphPercNaiveLE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64) <: DoubleGraph
#
# TODO
# """
function GraphPercNaiveLE(N::Integer, P::Integer, M::Integer, γ::Float64, β::Float64)
    ξ, ξv = PercNaive.gen_ξ(N, P)
    GraphLocEntr(N, M, γ, β, GraphPercNaive, ξ, ξv)
end

function GraphPercNaiveLE(X::GraphPercNaive, M::Integer, γ::Float64, β::Float64)
    GraphLocEntr(X.N, M, γ, β, GraphPercNaive, X.ξ, X.ξv)
end

end # module
