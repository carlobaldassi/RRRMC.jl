module QAliases

using ..QT
using ..Empty
using ..IsingSK
using ..EA

export GraphQ0T, GraphQIsingT, GraphQEAT

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
The graph is fully-connected, the interactions are random (\$J ∈ {-1,1}\$),
there are no external longitudinal fields.

See also [`Qenergy`](@ref).
"""
GraphQIsingT(Nk::Integer, M::Integer, Γ::Float64, β::Float64) = GraphQuant(Nk, M, Γ, β, GraphIsingSK, IsingSK.gen_J(Nk))




typealias GraphQEAT{fourK,twoD} GraphQuant{fourK,GraphEAContSimple{twoD}}

# """
#     GraphQIsingT(N::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph
#
# A `DoubleGraph` which implements a quantum Ising spin model in a transverse magnetic field,
# using the Suzuki-Trotter transformation.
# `N` is the number of spins, `M` the number of Suzuki-Trotter replicas, `Γ` the transverse
# field, `β` the inverse temperature.
# The graph is fully-connected, the interactions are random (\$J ∈ {-1,1}\$),
# there are no external longitudinal fields.
#
# See also [`Qenergy`](@ref).
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

end # module
