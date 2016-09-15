module LE

using ExtractMacro
using ..Interface
using ..Common

export GraphLE, GraphLocEntr, LEenergies

import ..Interface: energy, delta_energy, neighbors, allΔE,
                    update_cache!, delta_energy_residual

import Base: start, next, done, length, eltype

type GraphLE{M,γT} <: DiscrGraph{Float64}
    N::Int
    Nk::Int
    cache::LocalFields{Int}

    function GraphLE(N::Integer)
        isa(M, Int) || throw(ArgumentError("invalid parameter M, expected Int, given: $(typeof(M))"))
        M > 2 || throw(ArgumentError("M must be greater than 2, given: $M"))
        isa(γT, Float64) || throw(ArgumentError("invalid parameter γT, expected Float64, given: $(typeof(γT))"))
        N % (M + 1) == 0 || throw(ArgumentError("N must be divisible by M+1, given: N=$N M=$M"))
        Nk = N ÷ (M+1)

        cache = LocalFields{Int}(N)
        return new(N, Nk, cache)
    end
end

@doc """
    GraphLE{M,γT}(N::Integer) <: DiscrGraph

    TODO
""" -> GraphLE{M,γT}(N::Integer)

GraphLE{M,oldγ}(X::GraphLE{M,oldγ}, newγ::Float64) = GraphLE{M,newγ}(X.N)

# @generated function ΔElist{M,γT}(::Type{GraphLE{M,γT}})
#     Expr(:tuple, ntuple(d->fk(2*(d - 1 - (M-1) >>> 0x1) - iseven(M), γT), M)...)
# end
# lstind(μ̄::Int, M::Int) = (μ̄ + M-1) >>> 0x1 + 1
#
# function getk{M,γT}(X::GraphLE{M,γT}, μ̄::Int)
#     @inbounds k = ΔElist(GraphLE{M,γT})[lstind(μ̄, M)]
#     return k
# end

function energy{M,γT}(X::GraphLE{M,γT}, C::Config)
    # @assert X.N == C.N
    @extract X : Nk cache
    @extract cache : lfields lfields_last
    @extract C : s

    n = 0
    j = 0
    for i = 1:Nk
        j += 1
        # @assert j == 1 + (i-1) * (M+1)
        jc = j
        σc = 2s[jc] - 1

        μ = 0
        for k = 2:(M+1)
            j += 1
            # @assert j == k + (i-1) * (M+1)
            σj = 2s[j] - 1
            lfields[j] = σc * σj
            μ += σj
        end
        f = σc * μ
        lfields[jc] = f
        n -= f
    end
    cache.move_last = 0
    fill!(lfields_last, 0.0)
    return n * γT
end

function kinterval(move::Integer, M::Integer)
    j0 = move - ((move-1) % (M+1))
    j1 = j0 + M
    return j0:j1
end

function update_cache!{M,γT}(X::GraphLE{M,γT}, C::Config, move::Int)
    # @assert X.N == C.N
    # @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : Nk cache

    k = mod1(move, M+1)
    if k ≠ 1
        @inbounds σx = 2s[move] - 1
        i = (move-1) ÷ (M+1) + 1
        jc = 1 + (i-1) * (M+1)
    else
        Ux = kinterval(move, M)
    end

    @extract cache : lfields lfields_last move_last
    if move_last == move
        @inbounds begin
            #for y in neighbors(X, move)
            if k ≠ 1
                lfields[jc], lfields_last[jc] = lfields_last[jc], lfields[jc]
                lfields[move], lfields_last[move] = lfields_last[move], lfields[move]
            else
                for y in Ux
                    lfields[y], lfields_last[y] = lfields_last[y], lfields[y]
                end
            end
            #lfields[move] = -lfields[move]
            #lfields_last[move] = -lfields_last[move]
        end
        return
    end

    @inbounds begin
        if k ≠ 1
            #jc = 1 + (i-1) * (M+1)
            σc = 2s[jc] - 1
            lfc = lfields[jc]
            lfm = lfields[move]
            lfields_last[jc] = lfc
            lfields_last[move] = lfm
            lfields[jc] = lfc + 2 * (σc * σx)
            lfields[move] = -lfm
        else
            for y in Ux
                lfy = lfields[y]
                lfields_last[y] = lfy
                lfields[y] = -lfy
            end
        end
        #lfm = lfields[move]
        #lfields_last[move] = lfm
        #lfields[move] = -lfm
    end
    cache.move_last = move

    # lfields_bk = copy(lfields)
    # energy(X, C)
    # lfields_bk ≠ lfields && @show move hcat(lfields,lfields_bk) find(lfields .≠ lfields_bk)
    # @assert lfields_bk == lfields

    return
end

@inline function delta_energy{M,γT}(X::GraphLE{M,γT}, C::Config, move::Int)
    # @assert X.N == C.N
    # @assert 1 ≤ move ≤ C.N
    @extract X : cache
    @extract cache : lfields

    @inbounds Δ = 2γT * lfields[move]
    return Δ
end

immutable CavityRange
    j0::Int
    j1::Int
    jX::Int
    function CavityRange(j0::Integer, j1::Integer, jX::Integer)
        j0 ≤ jX ≤ j1 || throw(ArgumentError("invalid CavityRange parameters, expected j0≤jX≤j1, given: j0=$j0, j1=$j1, jX=$X"))
        return new(j0, j1, jX)
    end
end

start(crange::CavityRange) = crange.j0 + (crange.jX == crange.j0)
done(crange::CavityRange, j) = j > crange.j1
@inline function next(crange::CavityRange, j)
    @extract crange : j0 j1 jX
    # @assert j ≠ jX
    nj = j + 1
    nj += (nj == jX)
    return (j, nj)
end
length(crange::CavityRange) = crange.j1 - crange.j0
eltype(::Type{CavityRange}) = Int

@inline function neighbors{M}(X::GraphLE{M}, j::Int)
    j0 = j - ((j-1) % (M+1))
    j1 = j0 + M
    return CavityRange(j0, j1, j)
end

@generated function allΔE{M,γT}(::Type{GraphLE{M,γT}})
    iseven(M) ? Expr(:tuple, insert!([4*d*γT for d = 0:(M÷2)], 2, 2γT)...) :
                Expr(:tuple, [2*(2d-1)*γT for d = 1:((M+1)÷2)]...)
end

# Replicate an existsing graph

type GraphLocEntr{M,γT,G<:AbstractGraph} <: DoubleGraph{Float64}
    N::Int
    Nk::Int
    X0::GraphLE{M,γT}
    X1::Vector{G}
    C1::Vector{Config}
    function GraphLocEntr(N::Integer, g0::G, Gconstr, args...)
        X0 = GraphLE{M,γT}(N)
        Nk = X0.Nk
        X1 = Array{G}(M)
        X1[1] = g0
        for k = 2:M
            X1[k] = Gconstr(args...)
        end
        C1 = [Config(Nk, init=false) for k = 1:M]
        return new(N, Nk, X0, X1, C1)
    end
end

#  """
#      GraphLocEntr(...)
#
#  TODO
#  """
function GraphLocEntr(Nk::Integer, M::Integer, γ::Float64, β::Float64, Gconstr, args...)
    g0 = Gconstr(args...)
    G = typeof(g0)
    return GraphLocEntr{M,γ/β,G}(Nk * (M+1), g0, Gconstr, args...)
end

function update_cache!{M}(X::GraphLocEntr{M}, C::Config, move::Int)
    @extract X : X0 X1 C1
    k = mod1(move, M+1)
    i = (move - 1) ÷ (M+1) + 1

    k > 1 && spinflip!(X1[k-1], C1[k-1], i)

    update_cache!(X0, C, move)
end

function energy{M}(X::GraphLocEntr{M}, C::Config)
    # @assert X.N == C.N
    @extract X : Nk X0 X1 C1
    @extract C : s

    E = energy(X0, C)

    for k = 1:M
        s1 = C1[k].s
        for (i,j) = enumerate((k+1):(M+1):(k+1 + (M+1) * (Nk-1)))
            s1[i] = s[j]
        end
        E += energy(X1[k], C1[k])
    end

    return E
end

function LEenergies{M}(X::GraphLocEntr{M})
    @extract X : X1 C1

    Es = zeros(M)

    for k = 1:M
        Es[k] = energy(X1[k], C1[k])
    end

    return Es
end

function delta_energy_residual{M}(X::GraphLocEntr{M}, C::Config, move::Int)
    @extract X : X1 C1

    k = mod1(move, M+1)
    k == 1 && return 0.0

    i = (move - 1) ÷ (M+1) + 1

    return delta_energy(X1[k-1], C1[k-1], i)
end

function delta_energy(X::GraphLocEntr, C::Config, move::Int)
    return delta_energy(X.X0, C, move) +
           delta_energy_residual(X, C, move)
end

end
