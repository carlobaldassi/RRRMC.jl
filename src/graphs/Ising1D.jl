module Ising1D

using ExtractMacro
using ..Interface

export GraphIsing1D

import ..Interface: energy, delta_energy, neighbors, allΔE

type GraphIsing1D <: DiscrGraph{Int}
    N::Int
    J::BitVector
    sJ::Int
    function GraphIsing1D(N::Integer)
        @assert N > 2
        #J = bitrand(N)
        J = trues(N)
        sJ = sum(J)
        return new(N, J, sJ)
    end
end

@doc """
    GraphIsing1D(N::Integer) <: DiscrGraph

A simple 1-dimensional `DiscrGraph` type with `N` spins, antiferromagnetic interactions
(\$J=-1\$), no fields, and periodic boundary conditions.

Mostly useful for testing/debugging purposes.
""" -> GraphIsing1D(N::Integer)


function energy(X::GraphIsing1D, C::Config)
    @assert X.N == C.N
    @extract C : N s
    @extract X : J sJ

    #=n0 = 0
    for i = 1:N
        n0 += (2J[i] - 1) * (2s[i] - 1) * (2s[mod1(i+1,N)] - 1)
    end=#

    s1 = rol(s, 1)

    Js = J & s
    Js1 = J & s1
    ss1 = s & s1

    n1 = 8 * sum(Js & s1) - 4 * sum(Js) - 4 * sum(Js1) - 4 * sum(ss1) + 2 * sJ + 4 * sum(s) - N

    n1 += N - 2 * sum(s) # !!!

    #@assert n0 == n1
    return n1
end

function delta_energy(X::GraphIsing1D, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : J

    #=oldn = energy(X, C)
    s[move] $= 1
    newn = energy(X, C)
    s[move] $= 1
    Δn0 = newn - oldn=#

    Δn = 0
    if move == 1
        Δn -= (2J[N] - 1) * (2s[N] - 1) * (2s[1] - 1)
    else
        Δn -= (2J[move-1] - 1) * (2s[move-1] - 1) * (2s[move] - 1)
    end
    if move == N
        Δn -= (2J[N] - 1) * (2s[N] - 1) * (2s[1] - 1)
    else
        Δn -= (2J[move] - 1) * (2s[move] - 1) * (2s[move+1] - 1)
    end
    Δn *= 2

    Δn += 2 * (2 * s[move] - 1) # !!!

    #@assert Δn == Δn0 (Δn,Δn0)

    return Δn
end

#neighbors(X::GraphIsing1D, i::Int) = (i > 1 ? i-1 : X.N), (i < X.N ? i + 1 : 1)
neighbors(X::GraphIsing1D, i::Int) = (mod1(i-1, X.N), mod1(i+1, X.N))
allΔE(::Type{GraphIsing1D}) = (2,  6)



end

