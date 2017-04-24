# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

# This implements the functions needed by the waiting time method
# for the choice of the next move

module WaitingTimes

using ExtractMacro
using DataStructures
using ..Common
using ..Interface

export THeap, pick_next, update_heap!

τΔE(X::AbstractGraph, C::Config, β::Real, ΔE) = max(1.0, Float64(exp(β * ΔE)))
τ(X::AbstractGraph, C::Config, β::Real, i::Integer) = τΔE(X, C, β, delta_energy(X, C, i))
gen_wt(τ::Float64) = begin
    wt = -τ * log1p(-rand())
    #@assert wt ≥ 0 wt, τ, log1p(-rand())
    return wt
end

const THeap = MutableBinaryHeap{Float64, DataStructures.LessThan}

function THeap(X::AbstractGraph, C::Config, β::Real)
    N = getN(X)
    @assert C.N == N
    τs = [τ(X, C, β, i) for i = 1:N]
    theap = mutable_binary_minheap(Float64)
    for i = 1:N
        j = push!(theap, gen_wt(τs[i]))
        @assert j == i # implementation detail check
    end
    return theap
end

function pick_next(theap::THeap)
    # warning: depends on implementation details
    # see issue #217 in DataStructures.jl
    el = theap.nodes[1]
    t, i = el.value, el.handle
    return t, i
end

function update_heap!(theap::THeap, X::AbstractGraph, C::Config, β::Real, i::Integer, t::Float64)
    ΔE = delta_energy(X, C, i)
    spinflip!(X, C, i)
    τi = τΔE(X, C, β, -ΔE)
    ti = t + gen_wt(τi)
    update!(theap, i, ti)
    for j in neighbors(X, i)
        τj = τ(X, C, β, j)
        tj = t + gen_wt(τj)
        update!(theap, j, tj)
    end
    return ΔE
end

end # module
