module RRG

using ExtractMacro
using ..Interface
using ..Common
using ..DFloats

if isdefined(Main, :Documenter)
# this is silly but it's required for correct cross-linking in docstrings, apparently
using ...RRRMC
end

export GraphRRG, GraphRRGCont, GraphRRGContSimple

import ..Interface: energy, delta_energy, neighbors, allΔE,
                    delta_energy_residual, update_cache!, update_cache_residual!

import ..DFloats: MAXDIGITS
sentinel{ET}(::Type{ET}) = typemin(ET)

discr{ET}(::Type{ET}, x::Real) = convert(ET, round(x, MAXDIGITS))
discr(::Type{DFloat64}, x::Real) = x
discr{ET<:Integer}(::Type{ET}, x::Integer) = convert(ET, x)

# Bollobás pairing model generator for random regular graphs.
# Will fail when `K` is too large.
function gen_RRG(N::Integer, K::Integer)
    @assert K ≥ 1
    iseven(N * K) || throw(ArgumentError("N * K must be even, given N=$N, K=$K"))

    NK = N * K
    rp = Array(Int, NK)
    B = [falses(N) for i = 1:N]

    maxattmpts = 100_000
    ok = false
    for it = 1:maxattmpts
        resize!(rp, NK)
        rp[:] = 1:NK
        map(b->fill!(b, 0), B)
        l = NK
        again = false
        while !isempty(rp)
            #@show rp
            l = length(rp)
            j = rand(1:(l-1))
            rv1 = pop!(rp)
            rp[j], rp[end] = rp[end], rp[j]
            rv2 = pop!(rp)
            v1 = mod1(rv1, N)
            v2 = mod1(rv2, N)
            if v1 == v2 || B[v1][v2]
                #println("again!")
                again = true
                break
            end
            B[v1][v2] = 1
            B[v2][v1] = 1
        end
        again && continue
        ok = true
        break
    end
    ok || error("failed!")
    A = [find(b) for b in B]
    @assert all(a->length(a) == K, A)

    tA = NTuple{K,Int}[tuple(Ax...) for Ax in A]
    return tA
end

function gen_J{K}(f, ET::Type, N::Integer, A::Vector{NTuple{K,Int}})
    m = sentinel(ET)
    J = Vector{ET}[zeros(ET,K) for i = 1:N]
    map(Jx->fill!(Jx, m), J)
    for x = 1:N
        Jx = J[x]
        Ax = A[x]
        for k = 1:length(Ax)
            y = Ax[k]
            if x < y
                Jxy = f() #rand(vLEV)
                @assert Jx[k] == m
                Jx[k] = Jxy
                J[y][findfirst(J[y],m)] = Jxy
            else
                l = findfirst(A[y], x)
                Jxy = J[y][l]
                @assert Jxy ≠ m
                @assert J[x][k] == Jxy
            end
        end
    end
    @assert all(Jx->all(Jx .≠ m), J)
    tJ = NTuple{K,ET}[tuple(Jx...) for Jx in J]
    return tJ
end

function get_vLEV(LEV, ET::Type)
    isa(LEV, Tuple{Real,Vararg{Real}}) || throw(ArgumentError("invalid level spec, expected a Tuple of Reals, given: $LEV"))
    length(unique(LEV)) == length(LEV) || throw(ArgumentError("repeated levels in LEV: $LEV"))

    vLEV = Array(ET, length(LEV))
    try
        for i = 1:length(LEV)
            vLEV[i] = LEV[i]
        end
    catch
        throw(ArgumentError("incompatible energy type and level spec, conversion failed, given: ET=$ET LEV=$LEV"))
    end
    m = sentinel(ET)
    m ∈ vLEV && throw(ArgumentError("illegal level value: $m"))
    any(x->discr(ET, x) ≠ x, vLEV) && throw(ArgumentError("up to $MAXDIGITS decimal digits supported in levels, given: $LEV"))
    return vLEV
end

type GraphRRG{ET,LEV,K} <: DiscrGraph{ET}
    N::Int
    A::Vector{NTuple{K,Int}}
    J::Vector{NTuple{K,ET}}
    uA::Vector{Vector{Int}}
    cache::LocalFields{ET}
    function GraphRRG(A::Vector{NTuple{K,Int}}, J::Vector{NTuple{K,ET}})
        isa(K, Integer) || throw(ArgumentError("K must be integer, given a: $(typeof(K))"))

        N = length(A)
        #all(a->length(a) == K, A) || throw(ArgumentError("invalid A inner length, expected $K, given: $(unique(map(length,A)))"))

        vLEV = get_vLEV(LEV, ET)
        all(Jx->all(Jxy->Jxy ∈ vLEV, Jx), J) || throw(ArgumentError("the given J is incompatible with levels $LEV"))
        length(J) == N || throw(ArgumentError("incompatible lengths of A and J: $(length(A)), $(length(J))"))
        #all(j->length(j) == K, J) || throw(ArgumentError("invalid J inner length, expected $K, given: $(unique(map(length,J)))"))

        uA = [collect(A[x][find(J[x])]) for x = 1:N]

        cache = LocalFields{ET}(N)

        return new(N, A, J, uA, cache)
    end
end

"""
    GraphRRG(N::Integer, K::Integer, LEV = (-1,1)) <: DiscrGraph

A `DiscGraph` implementing a random regular graph with `N` spins and connectivity `K`.
*Note*: `N*K` must be even. Also, the graph generator uses the pairing model method by Bollobás,
with a cutoff on the number of restarts, and thus it may occasionally fail if `K` is large.
The interactions are extracted at random from `LEV`, which must be a `Tuple` of `Real`s.
No external fields.
"""
function GraphRRG{ET<:Real}(N::Integer, K::Integer, LEV::Tuple{ET,Vararg{ET}})
    A = gen_RRG(N, K)
    vLEV = get_vLEV(LEV, ET)
    J = gen_J(ET, N, A) do
        rand(vLEV)
    end
    return GraphRRG{ET,LEV,K}(A, J)
end
GraphRRG(N::Integer, K::Integer, LEV::Tuple{Real,Vararg{Real}}) = GraphRRG(N, K, promote(LEV...))
GraphRRG(N::Integer, K::Integer) = GraphRRG(N, K, (-1,1))

GraphRRG(N::Integer, K::Integer, LEV::Tuple{Float64,Vararg{Float64}}) = GraphRRG(N, K, map(DFloat64, LEV))

function energy{ET}(X::GraphRRG{ET}, C::Config)
    @assert X.N == C.N
    @extract C : s
    @extract X : A J cache
    @extract cache : lfields lfields_last
    n = zero(ET)
    for x = 1:length(A)
        Jx = J[x]
        σx = 2 * s[x] - 1
        Ax = A[x]
        lf = zero(ET)
        for k = 1:length(Ax)
            y = Ax[k]
            σy = 2 * s[y] - 1
            Jxy = Jx[k]

            lf -= Jxy * σx * σy
        end
        n += lf
        lfields[x] = discr(ET, 2lf)
    end
    n /= 2
    cache.move_last = 0
    fill!(lfields_last, zero(ET))
    return discr(ET, n)
end

function update_cache!{ET}(X::GraphRRG{ET}, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : A J cache

    @extract cache : lfields lfields_last move_last
    if move_last == move
        @inbounds begin
            Ax = A[move]
            for k = 1:length(Ax)
                y = Ax[k]
                lfields[y], lfields_last[y] = lfields_last[y], lfields[y]
            end
            lfields[move] = -lfields[move]
            lfields_last[move] = -lfields_last[move]
        end
        return
    end

    @inbounds begin
        Jx = J[move]
        sx = s[move]
        Ax = A[move]
        for k = 1:length(Ax)
            y = Ax[k]
            σxy = 1 - 2 * (sx $ s[y])
            Jxy = Jx[k]
            lfy = lfields[y]
            lfields_last[y] = lfy
            lfields[y] = discr(ET, lfy - 4 * σxy * Jxy)
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

function delta_energy{ET}(X::GraphRRG{ET}, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    #@extract C : s
    @extract X : cache
    @extract cache : lfields

    @inbounds Δ = -lfields[move]
    return Δ

    # @inbounds begin
    #     Δ = zero(ET)
    #     Jx = J[move]
    #     σx = 2 * s[move] - 1
    #     Ax = A[move]
    #     for k = 1:length(Ax)
    #         y = Ax[k]
    #         σy = 2 * s[y] - 1
    #         Jxy = Jx[k]
    #         Δ += 2 * Jxy * σx * σy
    #     end
    # end
    # return discr(ET, Δ)
end

neighbors(X::GraphRRG, i::Int) = return X.uA[i]
@generated _allΔE{K}(::Type{GraphRRG{Int,(-1,1),K}}) =
    iseven(K) ? Expr(:tuple, ntuple(d->4*(d-1), K÷2+1)...) :
                Expr(:tuple, ntuple(d->2*(2d-1), (K+1)÷2)...)


@generated function allΔE{ET,LEV,K}(::Type{GraphRRG{ET,LEV,K}})
    list = Set{ET}()
    L = length(LEV)
    es = Set{ET}(zero(ET))
    for n = 1:K
        newes = Set{ET}()
        for n in es, l in LEV
            push!(newes, discr(ET, n + l))
            push!(newes, discr(ET, n - l))
        end
        es = newes
    end
    deltas = sort!(unique(map!(x->2 * abs(x), collect(es))))
    return Expr(:tuple, deltas...)
end

## continous weights

type GraphRRGCont{ET,LEV,K} <: DoubleGraph{Float64}
    N::Int
    X0::GraphRRG{ET,LEV,K}
    A::Vector{NTuple{K,Int}}
    rJ::Vector{NTuple{K,Float64}}
    cache::LocalFields{Float64}
    function GraphRRGCont(N::Integer)
        A = gen_RRG(N, K)
        cJ = gen_J(Float64, N, A) do
            randn()
        end

        dJ = Array(NTuple{K,ET}, N)
        rJ = Array(NTuple{K,Float64}, N)
        for (x, cJx) in enumerate(cJ)
            dJ[x], rJ[x] = discretize(cJx, LEV)
        end

        X0 = GraphRRG{ET,LEV,K}(A, dJ)
        cache = LocalFields{Float64}(N)
        return new(N, X0, A, rJ, cache)
    end
end

"""
    GraphRRGCont(N::Integer, K::Integer, LEV) <: DoubleGraph{Float64}

A `DoubleGraph` implementing a random regular graph with `N` spins and connectivity `K`.
*Note*: `N*K` must be even. Also, the graph generator uses the pairing model method by Bollobás,
with a cutoff on the number of restarts, and thus it may occasionally fail if `K` is large.
The interactions are extracted from a normal distribution with unit variance, and are then
discretized using the values in `LEV`, which must be a `Tuple` of `Real`s. No external fields.

Same as [`GraphRRGContSimple`](@ref), but it can be used with [`rrrMC`](@ref).
"""
GraphRRGCont{ET<:Real}(N::Integer, K::Integer, LEV::Tuple{ET,Vararg{ET}}) = GraphRRGCont{ET,LEV,K}(N)
GraphRRGCont(N::Integer, K::Integer, LEV::Tuple{Real,Vararg{Real}}) = GraphRRGCont(N, K, promote(LEV...))

GraphRRGCont(N::Integer, K::Integer, LEV::Tuple{Float64,Vararg{Float64}}) = GraphRRGCont{DFloat64,map(DFloat64,LEV),K}(N)

function energy(X::GraphRRGCont, C::Config)
    @assert X.N == C.N
    @extract X : X0 A J=rJ cache
    @extract C : N s
    @extract cache : lfields lfields_last

    E0 = energy(X0, C)

    E1 = 0.0
    for x = 1:length(A)
        Jx = J[x]
        σx = 2 * s[x] - 1
        Ax = A[x]
        lf = 0.0
        for k = 1:length(Ax)
            y = Ax[k]
            σy = 2 * s[y] - 1
            Jxy = Jx[k]

            lf -= Jxy * σx * σy
        end
        E1 += lf
        lfields[x] = 2lf
    end
    E1 /= 2
    cache.move_last = 0
    fill!(lfields_last, 0.0)

    return convert(Float64, E0 + E1)
end

function update_cache!{ET}(X::GraphRRGCont{ET}, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N

    if X.cache.move_last ≠ X.X0.cache.move_last
        update_cache!(X.X0, C, move)
        update_cache_residual!(X, C, move)
        return
    end

    @extract C : N s
    @extract X : A rJ X0 cache
    @extract X0 : J0=J cache0=cache

    @extract cache : lfields lfields_last move_last
    @extract cache0 : lfields0=lfields lfields_last0=lfields_last move_last0=move_last
    @assert move_last0 == move_last
    if move_last == move
        @inbounds begin
            Ax = A[move]
            for k = 1:length(Ax)
                y = Ax[k]
                lfields[y], lfields_last[y] = lfields_last[y], lfields[y]
                lfields0[y], lfields_last0[y] = lfields_last0[y], lfields0[y]
            end
            lfields[move] = -lfields[move]
            lfields_last[move] = -lfields_last[move]
            lfields0[move] = -lfields0[move]
            lfields_last0[move] = -lfields_last0[move]
        end
        return
    end

    @inbounds begin
        Jx0 = J0[move]
        Jx = rJ[move]
        sx = s[move]
        Ax = A[move]
        for k = 1:length(Ax)
            y = Ax[k]
            σxy = 1 - 2 * (sx $ s[y])

            Jxy0 = Jx0[k]
            lfy0 = lfields0[y]
            lfields_last0[y] = lfy0
            lfields0[y] = discr(ET, lfy0 - 4 * σxy * Jxy0)

            Jxy = Jx[k]
            lfy = lfields[y]
            lfields_last[y] = lfy
            lfields[y] = lfy - 4 * σxy * Jxy
        end
        lfm0 = lfields0[move]
        lfields_last0[move] = lfm0
        lfields0[move] = -lfm0

        lfm = lfields[move]
        lfields_last[move] = lfm
        lfields[move] = -lfm
    end
    cache0.move_last = move
    cache.move_last = move

    #lfields_bk = copy(lfields)
    #energy(X, C)
    #@assert lfields_bk == lfields

    return
end

function update_cache_residual!(X::GraphRRGCont, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : A J=rJ cache

    @extract cache : lfields lfields_last move_last
    if move_last == move
        @inbounds begin
            Ax = A[move]
            for k = 1:length(Ax)
                y = Ax[k]
                lfields[y], lfields_last[y] = lfields_last[y], lfields[y]
            end
            lfields[move] = -lfields[move]
            lfields_last[move] = -lfields_last[move]
        end
        return
    end

    @inbounds begin
        Jx = J[move]
        sx = s[move]
        Ax = A[move]
        for k = 1:length(Ax)
            y = Ax[k]
            σxy = 1 - 2 * (sx $ s[y])
            Jxy = Jx[k]
            lfy = lfields[y]
            lfields_last[y] = lfy
            lfields[y] = lfy - 4 * σxy * Jxy
        end
        lfm = lfields[move]
        lfields_last[move] = lfm
        lfields[move] = -lfm
    end
    cache.move_last = move

    return
end

function delta_energy_residual(X::GraphRRGCont, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    #@extract C : s
    @extract X : cache
    @extract cache : lfields

    @inbounds Δ = -lfields[move]
    return Δ

    # @inbounds begin
    #     Δ = 0.0
    #     Jx = J[move]
    #     σx = 2 * s[move] - 1
    #     Ax = A[move]
    #     for k = 1:length(Ax)
    #         y = Ax[k]
    #         σy = 2 * s[y] - 1
    #         Jxy = Jx[k]
    #         Δ += 2 * Jxy * σx * σy
    #     end
    # end
    # return Δ
end

function delta_energy(X::GraphRRGCont, C::Config, move::Int)
    ΔE0 = delta_energy(X.X0, C, move)
    ΔE1 = delta_energy_residual(X, C, move)
    return convert(Float64, ΔE0 + ΔE1)
end


type GraphRRGContSimple{K} <: SimpleGraph{Float64}
    N::Int
    A::Vector{NTuple{K,Int}}
    J::Vector{NTuple{K,Float64}}
    cache::LocalFields{Float64}
    function GraphRRGContSimple(N::Integer)
        A = gen_RRG(N, K)
        J = gen_J(Float64, N, A) do
            randn()
        end

        cache = LocalFields{Float64}(N)
        return new(N, A, J, cache)
    end
end

"""
    GraphRRGContSimple(N::Integer, K::Integer) <: SimpleGraph{Flaot64}

A `SimpleGraph` implementing a random regular graph with `N` spins and connectivity `K`.
*Note*: `N*K` must be even. Also, the graph generator uses the pairing model method by Bollobás,
with a cutoff on the number of restarts, and thus it may occasionally fail if `K` is large.
The interactions are extracted from a normal distribution with unit variance.

Same as [`GraphRRGCont`](@ref), but it's more efficient when used with [`standardMC`](@ref).
"""
GraphRRGContSimple(N::Integer, K::Integer) = GraphRRGContSimple{K}(N)

function energy(X::GraphRRGContSimple, C::Config)
    @assert X.N == C.N
    @extract X : A J cache
    @extract C : N s
    @extract cache : lfields lfields_last

    E1 = 0.0
    for x = 1:length(A)
        Jx = J[x]
        σx = 2 * s[x] - 1
        Ax = A[x]
        lf = 0.0
        for k = 1:length(Ax)
            y = Ax[k]
            σy = 2 * s[y] - 1
            Jxy = Jx[k]

            lf -= Jxy * σx * σy
        end
        E1 += lf
        lfields[x] = 2lf
    end
    E1 /= 2
    cache.move_last = 0
    fill!(lfields_last, 0.0)

    return E1
end

function update_cache!(X::GraphRRGContSimple, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : A J cache

    @extract cache : lfields lfields_last move_last
    if move_last == move
        @inbounds begin
            Ax = A[move]
            for k = 1:length(Ax)
                y = Ax[k]
                lfields[y], lfields_last[y] = lfields_last[y], lfields[y]
            end
            lfields[move] = -lfields[move]
            lfields_last[move] = -lfields_last[move]
        end
        return
    end

    @inbounds begin
        Jx = J[move]
        sx = s[move]
        Ax = A[move]
        for k = 1:length(Ax)
            y = Ax[k]
            σxy = 1 - 2 * (sx $ s[y])
            Jxy = Jx[k]
            lfy = lfields[y]
            lfields_last[y] = lfy
            lfields[y] = lfy - 4 * σxy * Jxy
        end
        lfm = lfields[move]
        lfields_last[move] = lfm
        lfields[move] = -lfm
    end
    cache.move_last = move

    return
end

function delta_energy(X::GraphRRGContSimple, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    #@extract C : s
    @extract X : cache
    @extract cache : lfields

    @inbounds Δ = -lfields[move]
    return Δ

    # @inbounds begin
    #     Δ = 0.0
    #     Jx = J[move]
    #     σx = 2 * s[move] - 1
    #     Ax = A[move]
    #     for k = 1:length(Ax)
    #         y = Ax[k]
    #         σy = 2 * s[y] - 1
    #         Jxy = Jx[k]
    #         Δ += 2 * Jxy * σx * σy
    #     end
    # end
    # return Δ
end

neighbors(X::GraphRRGContSimple, i::Int) = return X.A[i]

end
