# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module PSpin3

using Random
using ExtractMacro
using ..Interface
using ..Common

export GraphPSpin3

import ..Interface: energy, delta_energy, neighbors, allΔE, update_cache!

struct GraphPSpin3{K} <: DiscrGraph{Int}
    N::Int
    #allA::Vector{Int} # aux contiguous memory
    A::Vector{Vector{Int}} # neighbors
    U::Vector{Vector{Int}} # unique neighbors
    #J::Vector{Vector{Int}} # currently unused (all 1)
    cache::LocalFields{Int}
    function GraphPSpin3(N::Integer, K::Integer)
        @assert K ≥ 1

        N % 3 == 0 || throw(ArgumentError("N must be divisible by 3, given: $N"))
        A = Vector{Int}[zeros(Int, 2*K) for i = 1:N]
        #allA = zeros(Int, 2 * K * N)
        #A = Array{Vector{Int}}(undef, N)
        #for i = 1:N
            #A[i] = pointer_to_array(pointer(allA, (i-1) * 2K + 1), 2K)
        #end

        for k = 1:K
            l = randperm(N)
            for i = 1:3:(N-2)
                v1, v2, v3 = l[i], l[i+1], l[i+2]
                A[v1][2*(k-1) + 1] = v2
                A[v1][2*(k-1) + 2] = v3
                A[v2][2*(k-1) + 1] = v1
                A[v2][2*(k-1) + 2] = v3
                A[v3][2*(k-1) + 1] = v1
                A[v3][2*(k-1) + 2] = v2
            end
        end

        @assert all(a->!any(x->x==0, a), A)

        U = Vector{Int}[unique(a) for a in A]

        cache = LocalFields{Int}(N)

        return new{K}(N, A, U, cache)
    end
end

@doc """
    GraphPSpin3(N::Integer, K::Integer) <: DiscrGraph{Int}

A `DiscrGraph` implementing a \$p\$-spin regular graph with \$p=3\$. `N` is the number of spins, and must
be divisible by \$3\$; `K` is the connectivity. All interactions are set to \$J=1\$.
""" -> GraphPSpin3(N, K)

function energy(X::GraphPSpin3{K}, C::Config) where {K}
    @assert X.N == C.N
    @extract C : s
    @extract X : A cache #J
    @extract cache : lfields lfields_last
    n = 0
    for x = 1:length(A)
        #Jx = J[x]
        σx = 2 * s[x] - 1
        Ax = A[x]
        p = 1
        lf = 0
        for k = 1:K
            y = Ax[p]
            z = Ax[p+1]
            σy = 2 * s[y] - 1
            σz = 2 * s[z] - 1
            #Jxyz = Jx[k]

            lf -= σy * σz
            p += 2
        end
        lf *= σx
        n += lf
        lfields[x] = 2lf
    end
    @assert n % 3 == 0
    n ÷= 3
    cache.move_last = 0
    fill!(lfields_last, 0)
    return n
end

function update_cache!(X::GraphPSpin3{K}, C::Config, move::Int) where {K}
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : A U cache

    @extract cache : lfields lfields_last move_last
    if move_last == move
        @inbounds begin
            Ux = U[move]
            for y in Ux
                lfields[y], lfields_last[y] = lfields_last[y], lfields[y]
            end
            lfields[move] = -lfields[move]
            lfields_last[move] = -lfields_last[move]
        end
        return
    end

    @inbounds begin
        Ux = U[move]
        for y in Ux
            lfields_last[y] = lfields[y]
        end
        #Jx = J[move]
        σx = 2 * s[move] - 1
        Ax = A[move]
        p = 1
        for k = 1:K
            y = Ax[p]
            z = Ax[p+1]
            σy = 2 * s[y] - 1
            σz = 2 * s[z] - 1
            #Jxyz = Jx[k]
            Δ = 4 * σx * σy * σz
            lfields[y] -= Δ
            lfields[z] -= Δ
            p += 2
        end
        lfm = lfields[move]
        lfields_last[move] = lfm
        lfields[move] = -lfm
    end
    cache.move_last = move

    #lfields_bk = copy(lfields)
    #energy(X, C)
    #@assert lfields_bk == lfields

    return
end

function delta_energy(X::GraphPSpin3, C::Config, move::Int) where {K}
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    #@extract C : s
    @extract X : cache # A #J
    @extract cache : lfields

    @inbounds Δ = -lfields[move]
    return Δ

    # @inbounds begin
    #     Δ = 0
    #     #Jx = J[move]
    #     σx = 2 * s[move] - 1
    #     Ax = A[move]
    #     p = 1
    #     for k = 1:K
    #         y = Ax[p]
    #         z = Ax[p+1]
    #         σy = 2 * s[y] - 1
    #         σz = 2 * s[z] - 1
    #         #Jxyz = Jx[k]
    #         Δ += σy * σz
    #         p += 2
    #     end
    #     Δ *= 2 * σx
    # end
    # return Δ
end

neighbors(X::GraphPSpin3, i::Int) = return X.U[i]
@generated allΔE(::Type{GraphPSpin3{K}}) where {K} =
    iseven(K) ? Expr(:tuple, ntuple(d->4*(d-1), K÷2+1)...) :
                Expr(:tuple, ntuple(d->2*(2d-1), (K+1)÷2)...)

end
