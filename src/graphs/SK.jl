# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module SK

using Random, LinearAlgebra
using ExtractMacro
using ..Interface
using ..Common
using ..QT

using ...RRRMC # this is silly but it's required for correct cross-linking in docstrings, apparently

export GraphSK, GraphSKNormal

import ..Interface: energy, delta_energy, update_cache!, neighbors

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

struct GraphSK <: SimpleGraph{Float64}
    N::Int
    sN::Float64
    J::Vector{BitVector}
    #tmps::BitVector
    cache::LocalFields{Int}
    function GraphSK(J::Vector{BitVector}; check::Bool = true)
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
        #tmps = BitArray(undef, N)
        cache = LocalFields{Int}(N)
        return new(N, √N, J, cache)
    end
end

"""
    GraphSK(N::Integer) <: SimpleGraph{Float64}

A `SimpleGraph` implementing a Sherrington-Kirkpatrick fully-connected Ising model
with `N` spins and random binary interactions (\$J ∈ \\{-1/√N,1/√N\\}\$) and no
external fields.

Same as [`GraphSKNormal`](@ref), but with binary interactions.
"""
GraphSK(N::Integer) = GraphSK(gen_J(N), check=false)

function energy(X::GraphSK, C::Config)
    @assert X.N == C.N
    @extract C : s
    @extract X : N sN J cache
    @extract cache : lfields lfields_last
    tmps = BitArray(undef, N)
    n = -2 * sum(s)
    for i = 1:N
        Ji = J[i]
        sc = sum(map!(⊻, tmps, Ji, s))
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

    return n / sN
end

function update_cache!(X::GraphSK, C::Config, move::Int)
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
            Jσij = si ⊻ s[j] ⊻ Ji[j]
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

function delta_energy(X::GraphSK, C::Config, move::Int)
    @extract X : sN cache
    @extract cache : lfields

    @inbounds Δ = lfields[move] / sN
    return Δ

    # @extract X : N J tmps
    # @extract C : s
    # @assert N == length(s)
    # @assert 1 ≤ move ≤ N

    # Ji = J[move]
    # si = s[move]
    # sc = sum(map!(⊻, tmps, Ji, s)) - si
    # @assert Δ == 2 * (2si-1) * (N-1 - 2sc)
    # return Δ
end

# function check_delta(X::GraphSK, C::Config, move::Int)
#     @extract C : s
#     delta = delta_energy(X, s, move)
#     e0 = energy(X, s)
#     s[move] ⊻= 1
#     e1 = energy(X, s)
#     s[move] ⊻= 1
#
#     (e1-e0) == delta || (@show e1,e0,delta,e1-e0; error())
# end

neighbors(X::GraphSK, i::Int) = AllButOne(X.N, i)

#####


function gen_J_gauss(N::Integer)
    J = Vec[rmul!(randn(N), 1/√N) for i = 1:N]
    for i = 1:N
        J[i][i] = 0
        for j = (i+1):N
            J[j][i] = J[i][j]
        end
    end
    return J
end

struct GraphSKNormal <: SimpleGraph{Float64}
    N::Int
    J::Vec2
    cache::LocalFields{Float64}
    function GraphSKNormal(J::Vec2; check::Bool = true)
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
        cache = LocalFields{Float64}(N)
        return new(N, J, cache)
    end
end

"""
    GraphSK(N::Integer) <: SimpleGraph{Float64}

A `SimpleGraph` implementing a Sherrington-Kirkpatrick fully-connected Ising model
with `N` spins and random interactions extracted from a normal distribution with
zero mean and \$1/N\$ variance, and no external fields.

Same as [`GraphSK`](@ref), but with Gaussian interactions.
"""
GraphSKNormal(N::Integer) = GraphSKNormal(gen_J_gauss(N), check=false)

function energy(X::GraphSKNormal, C::Config)
    @assert X.N == C.N
    @extract C : s
    @extract X : N J cache
    @extract cache : lfields lfields_last

    n = 0.0
    for i = 1:N
        Ji = J[i]
        si = s[i]
        lf = 0.0
        for j = 1:N
            # j == i && continue # avoid branching
            #n -= (2s[j] - 1) * (2s[i] - 1) * Ji[j]
            lf += (1 - 2 * (si ⊻ s[j])) * Ji[j]
        end
        lfields[i] = 2lf
        n -= lf
    end
    n /= 2

    cache.move_last = 0
    fill!(lfields_last, 0)

    return n
end

function update_cache!(X::GraphSKNormal, C::Config, move::Int)
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
            Jσij = (1 - 2 * (si ⊻ s[j])) * Ji[j]
            lfj = lfields[j]
            lfields_last[j] = lfj
            lfields[j] = lfj + 4 * Jσij
        end
        lfields_last[move] = lfm
        lfields[move] = -lfm
    end
    cache.move_last = move

    # lfields_bk = copy(lfields)
    # lfields_last_bk = copy(lfields_last)
    # energy(X, C)
    # @assert maximum(abs.(lfields_bk - lfields)) < 1e-10
    # copy!(lfields_last, lfields_last_bk)
    # cache.move_last = move

    return
end

function delta_energy(X::GraphSKNormal, C::Config, move::Int)
    @extract X : cache
    @extract cache : lfields

    @inbounds Δ = lfields[move]
    return Δ
end

function check_delta(X::GraphSKNormal, C::Config, move::Int)
    @extract C : s
    delta = delta_energy(X, C, move)
    s[move] = s[move] ⊻ 1
    e1 = energy(X, C)
    s[move] = s[move] ⊻ 1
    e0 = energy(X, C)

    abs((e1-e0) - delta) < 1e-10 || (@show e1,e0,delta,e1-e0; error())
end

neighbors(X::GraphSKNormal, i::Int) = AllButOne(X.N, i)

end
