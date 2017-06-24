# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module KSAT

using ExtractMacro
using Compat
using ..Interface
using ..Common

if isdefined(Main, :Documenter)
# this is silly but it's required for correct cross-linking in docstrings, apparently
using ...RRRMC
end

export GraphKSAT

import ..Interface: energy, delta_energy, neighbors, allΔE,
                    delta_energy_residual, update_cache!, update_cache_residual!

function choose(N::Int, K::Int)
    out = fill(typemax(Int), K)
    @inbounds for k = 1:K
        out[k] = rand(1:(N-k+1))
        for l = 1:(k-1)
            if out[l] ≤ out[k]
                out[k] += 1
            end
        end
        for l = 1:(k-1)
            if out[l] > out[k]
                x = out[k]
                for j = k:-1:l+1
                    out[j] = out[j-1]
                end
                out[l] = x
                break
            end
        end
    end
    #@assert issorted(out)
    #@assert all(diff(out) .≠ 0)
    return out
end

function gen_randomKSAT(N::Integer, K::Integer, α::Real)
    N > 0 || throw(ArgumentError("N must be positive: $N"))
    K > 0 || throw(ArgumentError("K must be positive: $K"))
    α ≥ 0 || throw(ArgumentError("α must be non-negative: $α"))
    N ≥ K || throw(ArgumentError("N must not be less than K: $N < $K"))

    M = round(Int, α * N)
    A = Array{IVec}(M)
    J = [BitArray(K) for a = 1:M]
    for a = 1:M
        A[a] = choose(N, K)
        rand!(J[a])
    end
    return A, J
end

type ClauseCache
    M::Int
    K::Int
    S::IVec
    I::IVec2
    ClauseCache(M::Integer, K::Integer) = new(M, K, zeros(Int, M), IVec[zeros(Int, K) for a = 1:M])
end

function clear!(cache::ClauseCache)
    @extract cache : S I
    fill!(S, 0)
    for Ia in I
        fill!(Ia, 0)
    end
    return cache
end

type GraphKSAT <: DiscrGraph{Int}
    N::Int
    M::Int
    K::Int
    A::IVec2
    J::Vector{BitVector}
    T::IVec2
    neighb::IVec2
    max_conn::Int
    cache::ClauseCache
    function GraphKSAT(N::Integer, K::Integer, A::IVec2, J::Vector{BitVector})
        M = length(A)
        length(J) == M || throw(ArgumentError("Incompatible lengths of A and J: $M vs $(length(J))"))

        T = [Int[] for i = 1:N]
        for a = 1:M
            for i in A[a]
                push!(T[i], a)
            end
        end

        neighb = [Int[] for i = 1:N]
        for i in 1:N
            Ti = T[i]
            for a in Ti, j in A[a]
                if j ≠ i && j ∉ neighb[i]
                    push!(neighb[i], j)
                end
            end
        end

        # TODO: more input consistency checks?
        max_conn = maximum(map(length, T))
        cache = ClauseCache(M, K)
        return new(N, M, K, A, J, T, neighb, max_conn, cache)
    end
end

"""
  GraphKSAT(N::Integer, α::Real, K::Integer)

A `DiscGraph` implementing a random `K`-SAT graph with `N` spins and `αN` clauses.
"""
function GraphKSAT(N::Integer, K::Integer, α::Real)
    A, J = gen_randomKSAT(N, K, α)
    return GraphKSAT(N, K, A, J)
end

function export_cnf(X::GraphKSAT, filename::AbstractString)
    @extract X : N M A J
    open(filename, "w") do f
        println(f, "p cnf $N $M")
        for a = 1:M
            for (j,i) in zip(J[a],A[a])
                print(f, (2j-1) * i, " ")
            end
            println(f, "0")
        end
    end
end

function export_cnf(X::GraphKSAT, filename::AbstractString, decimate::Vector{Int})
    @extract X : N M A=deepcopy(A) J=deepcopy(J) T=deepcopy(T)

    j = 1
    while j ≤ length(decimate)
        v = decimate[j]
        s = v > 0
        i = abs(v)
        for a in T[i]
            isempty(A[a]) && continue
            k = findfirst(A[a], i)
            @assert k > 0
            if J[a][k] == s
                empty!(A[a])
            else
                length(A[a]) > 1 || error("contradiction")
                deleteat!(A[a], k)
                deleteat!(J[a], k)
                if length(A[a]) == 1
                    newv = A[a][1] * (2J[a][1]-1)
                    -newv ∈ decimate && (error("contradiction"))
                    newv ∉ decimate && push!(decimate, newv)
                    empty!(A[a])
                end
            end
        end
        empty!(T[i])
        j += 1
    end

    nM = sum(map(x->!isempty(x),A)) + length(decimate)

    open(filename, "w") do f
        println(f, "p cnf $N $nM")
        for a = 1:M
            isempty(A[a]) && continue
            for (j,i) in zip(J[a],A[a])
                print(f, (2j-1) * i, " ")
            end
            println(f, "0")
        end
        for v in decimate
            println(f, "$v 0")
        end
    end
end

function energy(X::GraphKSAT, C::Config)
    length(C) == X.N || throw(ArgumentError("different N: $(length(C)) $(X.N)"))
    @extract C : s
    @extract X : M K A J cache
    clear!(cache)
    @extract cache : S I

    n = 0
    @inbounds for a = 1:M
        Ja = J[a]
        Aa = A[a]
        sat = 0
        Ia = I[a]
        for k = 1:K
            i = Aa[k]
            si = s[i]
            Ji = Ja[k]
            Ji ⊻ si == 0 && (sat += 1; Ia[sat] = i)
        end
        S[a] = sat
        sat == 0 && (n += 1)
    end
    return n
end

function delta_energy(X::GraphKSAT, C::Config, move::Int)
    @extract C : s
    @extract X : T cache
    @extract cache : S I

    Δ = 0
    @inbounds for a in T[move]
        Sa = S[a]
        Ia = I[a]
        if Sa == 1 && (Ia[1] == move)
            Δ += 1
        elseif S[a] == 0
            Δ -= 1
        end
    end
    return Δ
end

function update_cache!(X::GraphKSAT, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : K T cache
    @extract cache : S I

    # Δ = 0
    @inbounds for a in T[move]
        Sa = S[a]
        Ia = I[a]
        if Sa == 0
            S[a] = 1
            Ia[1] = move
            # Δ -= 1
            continue
        end
        sat = false
        for k = 1:Sa
            Ia[k] ≠ move && continue
            for l = k:Sa-1
                Ia[l] = Ia[l+1]
            end
            Ia[Sa] = 0
            S[a] = Sa - 1
            # Sa == 1 && (Δ += 1)
            sat = true
            break
        end
        sat && continue
        S[a] = Sa + 1
        Ia[Sa + 1] = move
    end
    return
end

neighbors(X::GraphKSAT, i::Int) = return X.neighb[i]

# TODO: improve
allΔE(X::GraphKSAT) = tuple(0:X.max_conn...)

end
