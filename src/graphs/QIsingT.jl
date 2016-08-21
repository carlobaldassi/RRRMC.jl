module QIsingT

using ExtractMacro
using ..Interface
using ..Common
using ..QT

export GraphQIsingT, GraphIsingSK, Qenergy

import ..Interface: energy, delta_energy, update_cache!
import ..QT: Qenergy, MAXDIGITS

function gen_J(N::Integer)
    J = BitVector[bitrand(N) for i = 1:N]
    for i = 1:N
        J[i][i] = 0
        for j = (i+1):N
            J[j][i] = J[i][j]
        end
    end
    return J
end

type GraphIsingSK <: SimpleGraph{Int}
    N::Int
    J::Vector{BitVector}
    #tmps::BitVector
    cache::LocalFields{Int}
    function GraphIsingSK(J::Vector{BitVector}; check::Bool = true)
        N = length(J)
        if check
            all(Jx->length(Jx) == N, J) || throw(ArgumentError("invalid J inner length, expected $N, given: $(unique(map(length,J)))"))
            for i = 1:N
                J[i][i] == 0 || throw(ArgumentError("diagonal entries of J must be 0, found: J[$i][$i] = $(J[i][i])"))
                for j = (i+1):N
                    J[i][j] == J[j][i] || throw(ArgumentError("J must be symmetric, found: J[$i][$j] = $(J[i][j]), J[$j][$i] = $(J[j][i])"))
                end
            end
        end
        #tmps = BitArray(N)
        cache = LocalFields{Int}(N)
        return new(N, J, cache)
    end
end

function energy(X::GraphIsingSK, C::Config)
    @assert X.N == C.N
    @extract C : s
    @extract X : N J cache
    @extract cache : lfields lfields_last
    tmps = BitArray(N)
    n = -2 * sum(s)
    for i = 1:N
        Ji = J[i]
        sc = sum(map!($, tmps, Ji, s))
        si = s[i]
        lf = -(2si-1) * (N-1 - 2sc)
        lfields[i] = 2 * (-lf + 2si)
        n += lf
    end
    @assert n % 2 == 0
    n ÷= 2
    cache.move_last = 0
    fill!(lfields_last, 0)

    # altn = 0
    # for i = 1:N
    #     Ji = J[i]
    #     for j = 1:N
    #         j == i && continue
    #         altn -= (2Ji[j] - 1) * (2s[j] - 1) * (2s[i] - 1)
    #     end
    # end
    # altn ÷= 2
    # @assert n == altn

    return n
end

function update_cache!(X::GraphIsingSK, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : J cache

    @extract cache : lfields lfields_last move_last

    if move_last == move
        cache.lfields, cache.lfields_last = cache.lfields_last, cache.lfields
        return
    end

    @inbounds begin
        Ji = J[move]
        si = s[move]
        lfm = lfields[move]
        @simd for j = 1:N
            # note: we don't check move ≠ j to avoid branching
            Jσij = si $ s[j] $ Ji[j]
            lfj = lfields[j]
            lfields_last[j] = lfj
            lfields[j] = lfj + 8 * Jσij - 4
        end
        lfields_last[move] = lfm
        lfields[move] = -lfm
    end
    cache.move_last = move

    # lfields_bk = copy(lfields)
    # lfields_last_bk = copy(lfields_last)
    # energy(X, C)
    # @assert lfields_bk == lfields
    # copy!(lfields_last, lfields_last_bk)
    # cache.move_last = move

    return
end

function delta_energy(X::GraphIsingSK, C::Config, move::Int)
    @extract X : cache
    @extract cache : lfields

    @inbounds Δ = lfields[move]
    return Δ

    # @extract X : N J tmps
    # @extract C : s
    # @assert N == length(s)
    # @assert 1 ≤ move ≤ N

    # Ji = J[move]
    # si = s[move]
    # sc = sum(map!($, tmps, Ji, s)) - si
    # @assert Δ == 2 * (2si-1) * (N-1 - 2sc)
    # return Δ
end

function check_delta(X::GraphIsingSK, C::Config, move::Int)
    @extract C : s
    delta = delta_energy(X, s, move)
    e0 = energy(X, s)
    s[move] $= 1
    e1 = energy(X, s)
    s[move] $= 1

    (e1-e0) == delta || (@show e1,e0,delta,e1-e0; error())
end


type GraphQIsingT{fourK} <: DoubleGraph{Float64}
    N::Int
    M::Int
    Nk::Int
    X0::GraphQT{fourK}
    X1::Vector{GraphIsingSK}
    C1::Vector{Config}
    λ::Float64
    H0::Float64 # useless!!!
    β::Float64
    Γ::Float64
    function GraphQIsingT(N::Integer, M::Integer, H0::Real, β::Float64, Γ::Float64)
        X0 = GraphQT{fourK}(N, M)
        Nk = X0.Nk
        J = gen_J(Nk)
        X1 = [GraphIsingSK(J, check=false) for k = 1:M]
        C1 = [Config(Nk, init=false) for k = 1:M]
        λ = 1 / (M * √N)
        return new(N, M, Nk, X0, X1, C1, λ, H0, β, Γ)
    end
end

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
function GraphQIsingT(Nk::Integer, M::Integer, Γ::Float64, β::Float64)
    @assert Γ ≥ 0
    fourK = round(2/β * log(coth(β * Γ / M)), MAXDIGITS)
    H0 = Nk * M / 2β * log(sinh(2β * Γ / M) / 2) # useless!!!
    return GraphQIsingT{fourK}(Nk * M, M, H0, β, Γ)
end

function update_cache!(X::GraphQIsingT, C::Config, move::Int)
    @extract X : X1 Nk C1
    #@extract C : s
    k = (move - 1) ÷ Nk + 1
    i = mod1(move, Nk)

    spinflip!(X1[k], C1[k], i)

    #s1 = C1[k].s
    #copy!(s1, 1, s, (k-1)*Nk + 1, Nk)
    #@assert C1[k].s == s[((k-1)*Nk + 1):k*Nk]
end

function energy(X::GraphQIsingT, C::Config)
    @assert X.N == C.N
    @extract X : M Nk X0 X1 C1 λ H0
    @extract C : s

    E = energy(X0, C)

    for k = 1:M
        s1 = C1[k].s
        copy!(s1, 1, s, (k-1)*Nk + 1, Nk)
        E += energy(X1[k], C1[k]) * λ
    end

    return E - H0
end

function Qenergy(X::GraphQIsingT, C::Config)
    @assert X.N == C.N
    @extract X : M Nk X0 X1 C1 λ β Γ
    #@extract C : s

    E = -Γ * transverse_mag(X0, C, β)

    for k = 1:M
        #s1 = C1[k].s
        #copy!(s1, 1, s, (k-1)*Nk + 1, Nk)
        #@assert C1[k].s == s[((k-1)*Nk + 1):k*Nk]
        E += energy(X1[k], C1[k]) * λ / Nk
    end

    return E
end

function delta_energy_residual(X::GraphQIsingT, C::Config, move::Int)
    @extract X : Nk X1 C1 λ
    #@extract C : s

    k = (move - 1) ÷ Nk + 1
    #s1 = C1[k].s
    #copy!(s1, 1, s, (k-1)*Nk + 1, Nk)
    #@assert C1[k].s == s[((k-1)*Nk + 1):k*Nk]

    i = mod1(move, Nk)
    return delta_energy(X1[k], C1[k], i) * λ
end

function delta_energy(X::GraphQIsingT, C::Config, move::Int)
    return delta_energy(X.X0, C, move) +
           delta_energy_residual(X, C, move)
end

end
