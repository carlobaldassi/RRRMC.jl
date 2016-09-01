module QAliases

using ..QT
using ..IsingSK
using ..Empty

export GraphQIsingT, GraphQ0T

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
GraphQIsingT(Nk::Integer, M::Integer, Γ::Float64, β::Float64) = GraphQuant(Nk, M, Γ, β, GraphIsingSK, gen_J(Nk))

end # module
