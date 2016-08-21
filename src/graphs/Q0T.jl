module Q0T

using ExtractMacro
using ..Interface
using ..Common
using ..QT

export GraphQ0T, Qenergy

import ..Interface: energy, delta_energy
import ..QT: Qenergy

const MAXDIGITS = QT.MAXDIGITS

type GraphQ0T{fourK} <: DoubleGraph{Float64}
    N::Int
    M::Int
    Nk::Int
    X0::GraphQT{fourK}
    H0::Float64 # useless!!!
    β::Float64
    Γ::Float64
    function GraphQ0T(N::Integer, M::Integer, H0::Real, β::Float64, Γ::Float64)
        X0 = GraphQT{fourK}(N, M)
        Nk = X0.Nk
        return new(N, M, Nk, X0, H0, β, Γ)
    end
end

"""
    GraphQ0T(N::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph

A simple `DoubleGraph` which implements independent spins in a transverse magnetic field,
using the Suzuki-Trotter transformation.
`N` is the number of spins, `M` the number of Suzuki-Trotter replicas,
`Γ` the transverse field, `β` the inverse temperature.

Intended for testing/debugging purposes.

See also [`Qenergy`](@ref).
"""
function GraphQ0T(Nk::Integer, M::Integer, Γ::Float64, β::Float64)
    @assert Γ ≥ 0
    fourK = round(2/β * log(coth(β * Γ / M)), MAXDIGITS)
    H0 = Nk * M / 2β * log(sinh(2β * Γ / M) / 2) # useless!!!
    return GraphQ0T{fourK}(Nk * M, M, H0, β, Γ)
end

energy(X::GraphQ0T, C::Config) = energy(X.X0, C) - X.H0

function Qenergy(X::GraphQ0T, C::Config)
    @assert X.N == C.N
    @extract X : X0 β Γ

    return -Γ * transverse_mag(X0, C, β)
end

delta_energy_residual(X::GraphQ0T, C::Config, move::Int) = 0.0
delta_energy(X::GraphQ0T, C::Config, move::Int) = delta_energy(X.X0, C, move)

end # module
