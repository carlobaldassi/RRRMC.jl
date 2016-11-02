module EA

using ExtractMacro
using ..Interface
using ..Common
using ..DFloats

if isdefined(Main, :Documenter)
using ...RRRMC # this is silly but it's required for correct cross-linking in docstrings, apparently
end

export GraphEA, GraphEACont, GraphEAContSimple

import ..Interface: energy, delta_energy, neighbors, allΔE,
                    update_cache!, update_cache_residual!

import ..DFloats: MAXDIGITS
sentinel{ET}(::Type{ET}) = typemin(ET)

discr{ET}(::Type{ET}, x::Real) = convert(ET, round(x, MAXDIGITS))
discr(::Type{DFloat64}, x::Real) = x
discr{ET<:Integer}(::Type{ET}, x::Integer) = convert(ET, x)

function gen_EA(L::Integer, D::Integer)
    L ≥ 2 || throw(ArgumentError("L must be ≥ 2, given: $L"))
    D ≥ 1 || throw(ArgumentError("D must be ≥ 0, given: $D"))
    N = L^D
    A = [Int[] for x = 1:N]
    dims = ntuple(d->L, D)

    for cind in CartesianRange(dims)
        x = sub2ind(dims, cind.I...)
        for d in 1:D
            I1 = ntuple(k->(k==d ? mod1(cind.I[k] + 1, L) : cind.I[k]), D)
            y = sub2ind(dims, I1...)
            push!(A[x], y)
            push!(A[y], x)
        end
    end
    map!(sort!, A)
    tA = NTuple{2D,Int}[tuple(Ax...) for Ax in A]
    return tA
end

function gen_J{twoD}(f, ET::Type, N::Integer, A::Vector{NTuple{twoD,Int}})
    @assert all(issorted, A)
    m = sentinel(ET)
    J = Vector{ET}[zeros(ET,twoD) for i = 1:N]
    map(Jx->fill!(Jx, m), J)
    for x = 1:N
        Jx = J[x]
        Ax = A[x]
        for k = 1:length(Ax)
            y = Ax[k]
            if x < y
                Jxy = f()
                @assert Jx[k] == m
                Jx[k] = Jxy
                J[y][findfirst(J[y],m)] = Jxy
            # else # this check fails for L=2
            #     l = findfirst(A[y], x)
            #     Jxy = J[y][l]
            #     @assert Jxy ≠ m
            #     @assert J[x][k] == Jxy
            end
        end
    end
    @assert all(Jx->all(Jx .≠ m), J)
    tJ = NTuple{twoD,ET}[tuple(Jx...) for Jx in J]
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

type GraphEA{ET,LEV,twoD} <: DiscrGraph{ET}
    N::Int
    #L::Int
    A::Vector{NTuple{twoD,Int}}
    J::Vector{NTuple{twoD,ET}}
    uA::Vector{Vector{Int}}
    cache::LocalFields{ET}
    function GraphEA(A::Vector{NTuple{twoD,Int}}, J::Vector{NTuple{twoD,ET}})
        isa(twoD, Integer) || throw(ArgumentError("twoD must be integer, given a: $(typeof(twoD))"))
        iseven(twoD) || throw(ArgumentError("twoD must be even, given: $twoD"))
        D = twoD ÷ 2
        D ≥ 1 || throw(ArgumentError("D must be ≥ 0, given: $D"))

        N = length(A)
        N ≥ 2 || throw(ArgumentError("invalid A, expected length ≥ 2, found: $N"))
        all(a->length(a) == twoD, A) || throw(ArgumentError("invalid A inner length, expected $twoD, given: $(unique(map(length,A)))"))

        isL2 = (A[1][1] == A[1][2])
        all(a->(issorted(a) && unique(a) == collect(a[1:(1+isL2):end])), A) || throw(ArgumentError("invalid A, does not look like an EA graph"))

        uA = [collect(a[1:(1+isL2):end]) for a in A]

        vLEV = get_vLEV(LEV, ET)
        all(Jx->all(Jxy->Jxy ∈ vLEV, Jx), J) || throw(ArgumentError("the given J is incompatible with levels $LEV"))
        length(J) == N || throw(ArgumentError("incompatible lengths of A and J: $(length(A)), $(length(J))"))
        all(j->length(j) == twoD, J) || throw(ArgumentError("invalid J inner length, expected $twoD, given: $(unique(map(length,J)))"))

        cache = LocalFields{ET}(N)

        return new(N, A, J, uA, cache)
    end
end

"""
    GraphEA(L::Integer, D::Integer, LEV = (-1,1)) <: DiscrGraph

An Edwards-Anderson `DiscrGraph`: spins are arranged on a square lattice of size `L`
in `D` dimensions (i.e. there are \$L^D\$ total spins), with periodic boundary
conditions. The interactions are extracted at random from `LEV`, which must be
a `Tuple` of `Real`s. No external fields.
"""
function GraphEA{ET<:Real}(L::Integer, D::Integer, LEV::Tuple{ET,Vararg{ET}})
    A = gen_EA(L, D)
    vLEV = get_vLEV(LEV, ET)
    N = length(A)
    J = gen_J(ET, N, A) do
        rand(vLEV)
    end
    return GraphEA{ET,LEV,2D}(A, J)
end
GraphEA(L::Integer, D::Integer, LEV::Tuple{Real,Vararg{Real}}) = GraphEA(L, D, promote(LEV...))
GraphEA(L::Integer, D::Integer) = GraphEA(L, D, (-1,1))

GraphEA(L::Integer, D::Integer, LEV::Tuple{Float64,Vararg{Float64}}) = GraphEA(L, D, map(DFloat64, LEV))

function energy{ET}(X::GraphEA{ET}, C::Config)
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
    #@assert n % 2 == 0
    #n ÷= 2
    n /= 2
    cache.move_last = 0
    fill!(lfields_last, zero(ET))
    return discr(ET, n)
end

function update_cache!{ET}(X::GraphEA{ET}, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : A uA J cache

    @extract cache : lfields lfields_last move_last
    if move_last == move
        @inbounds begin
            Ux = uA[move]
            for y in Ux
                lfields[y], lfields_last[y] = lfields_last[y], lfields[y]
            end
            lfields[move] = -lfields[move]
            lfields_last[move] = -lfields_last[move]
        end
        return
    end

    @inbounds begin
        Ux = uA[move]
        for y in Ux
            lfields_last[y] = lfields[y]
        end
        Jx = J[move]
        sx = s[move]
        Ax = A[move]
        for k = 1:length(Ax)
            y = Ax[k]
            σxy = 1 - 2 * (sx $ s[y])
            Jxy = Jx[k]
            lfields[y] = discr(ET, lfields[y] - 4 * σxy * Jxy)
        end
        lfm = lfields[move]
        lfields_last[move] = lfm
        lfields[move] = -lfm
    end
    cache.move_last = move

    return
end

function delta_energy{ET}(X::GraphEA{ET}, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    #@extract C : s
    #@extract X : A J
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

neighbors(X::GraphEA, i::Int) = return X.uA[i]
@generated allΔE{twoD}(::Type{GraphEA{Int,(-1,1),twoD}}) = Expr(:tuple, ntuple(d1->(4 * (d1 - 1)), (twoD÷2)+1)...)

@generated function allΔE{ET,LEV,twoD}(::Type{GraphEA{ET,LEV,twoD}})
    list = Set{ET}()
    L = length(LEV)
    es = Set{ET}(zero(ET))
    for n = 1:twoD
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

type GraphEACont{ET,LEV,twoD} <: DoubleGraph{Float64}
    N::Int
    X0::GraphEA{ET,LEV,twoD}
    A::Vector{NTuple{twoD,Int}}
    rJ::Vector{NTuple{twoD,Float64}}
    uA::Vector{Vector{Int}}
    cache::LocalFields{Float64}
    function GraphEACont(L::Integer)
        isa(twoD, Integer) || throw(ArgumentError("twoD must be integer, given a: $(typeof(twoD))"))
        iseven(twoD) || throw(ArgumentError("twoD must be even, given: $twoD"))
        D = twoD ÷ 2
        D ≥ 1 || throw(ArgumentError("D must be ≥ 0, given: $D"))
        A = gen_EA(L, D)
        N = length(A)
        cJ = gen_J(Float64, N, A) do
            randn()
        end

        dJ = Array(NTuple{twoD,ET}, N)
        rJ = Array(NTuple{twoD,Float64}, N)
        for (x, cJx) in enumerate(cJ)
            dJ[x], rJ[x] = discretize(cJx, LEV)
        end

        X0 = GraphEA{ET,LEV,twoD}(A, dJ)
        cache = LocalFields{Float64}(N)
        return new(N, X0, A, rJ, X0.uA, cache)
    end
end

"""
    GraphEACont(L::Integer, D::Integer, LEV) <: DoubleGraph{Float64}

An Edwards-Anderson `DoubleGraph`: spins are arranged on a square lattice of size `L`
in `D` dimensions (i.e. there are \$L^D\$ total spins), with periodic boundary
conditions. The interactions are extracted at random from a normal distribution
with unit variance, and are then discretized using the values in `LEV`,
which must be a `Tuple` of `Real`s. No external fields.

Same as [`GraphEAContSimple`](@ref), but it can be used with [`rrrMC`](@ref).
"""
GraphEACont{ET<:Real}(L::Integer, D::Integer, LEV::Tuple{ET,Vararg{ET}}) = GraphEACont{ET,LEV,2D}(L)
GraphEACont(L::Integer, D::Integer, LEV::Tuple{Real,Vararg{Real}}) = GraphEACont(L, D, promote(LEV...))

GraphEACont(L::Integer, D::Integer, LEV::Tuple{Float64,Vararg{Float64}}) = GraphEACont{DFloat64,map(DFloat64,LEV),2D}(L)

function energy(X::GraphEACont, C::Config)
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

function update_cache!{ET}(X::GraphEACont{ET}, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N

    if X.cache.move_last ≠ X.X0.cache.move_last
        update_cache!(X.X0, C, move)
        update_cache_residual!(X, C, move)
        return
    end

    @extract C : N s
    @extract X : A uA rJ X0 cache
    @extract X0 : J0=J cache0=cache

    @extract cache : lfields lfields_last move_last
    @extract cache0 : lfields0=lfields lfields_last0=lfields_last move_last0=move_last
    @assert move_last0 == move_last
    if move_last == move
        @inbounds begin
            Ux = uA[move]
            for y in Ux
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
        Ux = uA[move]
        for y in Ux
            lfields_last[y] = lfields[y]
            lfields_last0[y] = lfields0[y]
        end
        Jx0 = J0[move]
        Jx = rJ[move]
        sx = s[move]
        Ax = A[move]
        for k = 1:length(Ax)
            y = Ax[k]
            σxy = 1 - 2 * (sx $ s[y])

            Jxy0 = Jx0[k]
            lfields0[y] = discr(ET, lfields0[y] - 4 * σxy * Jxy0)

            Jxy = Jx[k]
            lfields[y] -= 4 * σxy * Jxy
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

    return
end

function update_cache_residual!(X::GraphEACont, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : A uA J=rJ cache

    @extract cache : lfields lfields_last move_last
    if move_last == move
        @inbounds begin
            Ux = uA[move]
            for y in Ux
                lfields[y], lfields_last[y] = lfields_last[y], lfields[y]
            end
            lfields[move] = -lfields[move]
            lfields_last[move] = -lfields_last[move]
        end
        return
    end

    @inbounds begin
        Ux = uA[move]
        for y in Ux
            lfields_last[y] = lfields[y]
        end
        Jx = J[move]
        sx = s[move]
        Ax = A[move]
        for k = 1:length(Ax)
            y = Ax[k]
            σxy = 1 - 2 * (sx $ s[y])
            Jxy = Jx[k]
            lfields[y] -= 4 * σxy * Jxy
        end
        lfm = lfields[move]
        lfields_last[move] = lfm
        lfields[move] = -lfm
    end
    cache.move_last = move

    return
end

function delta_energy_residual(X::GraphEACont, C::Config, move::Int)
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

function delta_energy(X::GraphEACont, C::Config, move::Int)
    ΔE0 = delta_energy(X.X0, C, move)
    ΔE1 = delta_energy_residual(X, C, move)
    return convert(Float64, ΔE0 + ΔE1)
end


type GraphEAContSimple{twoD} <: SimpleGraph{Float64}
    N::Int
    A::Vector{NTuple{twoD,Int}}
    J::Vector{NTuple{twoD,Float64}}
    uA::Vector{Vector{Int}}
    cache::LocalFields{Float64}
    function GraphEAContSimple(L::Integer)
        isa(twoD, Integer) || throw(ArgumentError("twoD must be integer, given a: $(typeof(twoD))"))
        iseven(twoD) || throw(ArgumentError("twoD must be even, given: $twoD"))
        D = twoD ÷ 2
        D ≥ 1 || throw(ArgumentError("D must be ≥ 0, given: $D"))
        A = gen_EA(L, D)
        N = length(A)
        J = gen_J(Float64, N, A) do
            randn()
        end

        uA = [unique(a) for a in A] # needed for the case L=2

        cache = LocalFields{Float64}(N)
        return new(N, A, J, uA, cache)
    end
end

"""
    GraphEACont(L::Integer, D::Integer) <: SimpleGraph{Float64}

An Edwards-Anderson `SimpleGraph`: spins are arranged on a square lattice of size `L`
in `D` dimensions (i.e. there are \$L^D\$ total spins), with periodic boundary
conditions. The interactions are extracted at random from a normal distribution
with unit variance.

Same as [`GraphEACont`](@ref), but it's more efficient when used with [`standardMC`](@ref).
"""
GraphEAContSimple(L::Integer, D::Integer) = GraphEAContSimple{2D}(L)

function energy(X::GraphEAContSimple, C::Config)
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

function update_cache!(X::GraphEAContSimple, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : A uA J cache

    @extract cache : lfields lfields_last move_last
    if move_last == move
        @inbounds begin
            Ux = uA[move]
            for y in Ux
                lfields[y], lfields_last[y] = lfields_last[y], lfields[y]
            end
            lfields[move] = -lfields[move]
            lfields_last[move] = -lfields_last[move]
        end
        return
    end

    @inbounds begin
        Ux = uA[move]
        for y in Ux
            lfields_last[y] = lfields[y]
        end
        Jx = J[move]
        sx = s[move]
        Ax = A[move]
        for k = 1:length(Ax)
            y = Ax[k]
            σxy = 1 - 2 * (sx $ s[y])
            Jxy = Jx[k]
            lfields[y] -= 4 * σxy * Jxy
        end
        lfm = lfields[move]
        lfields_last[move] = lfm
        lfields[move] = -lfm
    end
    cache.move_last = move

    return
end

function delta_energy(X::GraphEAContSimple, C::Config, move::Int)
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

neighbors(X::GraphEAContSimple, i::Int) = return X.uA[i]

end
