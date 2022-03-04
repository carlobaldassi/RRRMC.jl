# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

# This implements an efficient sampler from a discrete probability distribution
# which changes in time.
# See Wong and Easton 1980: http://dx.doi.org/10.1137/0209009

module DynamicSamplers

using Random
using ExtractMacro
using ..Common

import Base: getindex, setindex!, @propagate_inbounds, copy!, copy
import Random: rand

export DynamicSampler, refresh!, getel

mutable struct DynamicSampler
    v::Vec
    ps::Vec
    z::Float64
    N::Int
    levs::Int
    tinds::Vector{Int}
    tpos::Vector{Int}
    trefresh::Int
    function DynamicSampler(N::Int)
        levs = ceil(Int, log2(N))
        N2 = 2^levs
        v2 = zeros(N2)
        ps = zeros(N2 - 1)
        tinds, tpos = buildtable(levs)
        return new(v2, ps, 0.0, N, levs, tinds, tpos, 0)
    end
    function DynamicSampler(v)
        isempty(v) && throw(ArgumentError("input must be non-empty"))
        N = length(v)
        levs = ceil(Int, log2(N))
        N2 = 2^levs
        v2 = zeros(N2)
        v2[1:N] .= v
        n = findfirst(v2 .< 0)
        n ≡ nothing || throw(ArgumentError("negative element given: v[$n] = $(v2[n])"))
        ps = zeros(N2 - 1)
        tinds, tpos = buildtable(levs)

        dynsmp = new(v2, ps, 0.0, N, levs, tinds, tpos, 0)
        refresh!(dynsmp)

        return dynsmp
    end
end

function buildtable(levs::Int)
    N2 = 2^levs

    tinds = Array{Int}(undef, N2+1)
    tpos = empty!(Array{Int}(undef, N2 * levs))

    j0 = 1
    for i = 1:N2
        tinds[i] = j0
        j1 = j0
        k = 0
        u = 1 << (levs-1)
        off = 1
        @inbounds for lev = 1:levs
            if (i-1) & u == 0
                push!(tpos, off + k)
                j1 += 1
                k *= 2
            else
                k = 2k + 1
            end
            u >>>= 1
            off *= 2
        end
        j0 = j1
    end
    tinds[N2+1] = j0
    return tinds, tpos
end

function refresh!(dynsmp::DynamicSampler)
    @extract dynsmp : v ps N levs tinds tpos

    dynsmp.z = sum(v)
    fill!(ps, 0.0)
    @inbounds for i = 1:N
        x = v[i]
        @simd for j = tinds[i]:(tinds[i+1]-1)
            k = tpos[j]
            ps[k] += v[i]
        end
    end
    dynsmp.trefresh = 0
    return dynsmp
end

function copy!(dest::DynamicSampler, src::DynamicSampler)
    @extract src : v ps z N levs
    dest.N == N || throw(ArgumentError("dest N=$(dest.N), src N=$N"))
    dest.levs == levs || throw(ArgumentError("dest levs=$(dest.levs), src levs=$levs"))

    copy!(dest.v, v)
    copy!(dest.ps, ps)
    dest.z = src.z
    return dest
end

copy(dynsmp::DynamicSampler) = copy!(DynamicSampler(dynsmp.N), dynsmp)

# for checking purposes
function getel_naive(v::Vec, z::Float64, x::Float64)
    @assert 0 ≤ x < 1
    #@extract dynsmp : v z N
    N = length(v)
    x *= z
    i = 0
    c = 0.0
    @inbounds for i = 1:N
        c += v[i]
        c ≥ x && break
    end
    return i
end

getel_naive(dynsmp::DynamicSampler, x::Float64) = getel_naive(dynsmp.v, dynsmp.z, x)

function getel(dynsmp::DynamicSampler, x::Float64)
    #@assert 0 ≤ x < 1
    @extract dynsmp : v ps z N levs
    x *= z

    k = 0
    off = 1
    @inbounds for lev = 1:levs
        p = ps[off+k]
        k *= 2
        if x > p
            x -= p
            k += 1
        end
        off *= 2
    end
    if k ≥ N || v[k+1] == 0
        dynsmp.trefresh > 0 || throw(Error("Unrecoverable loss of precision detected in the dynamic sampler"))
        refresh!(dynsmp)
        return getel(dynsmp, x)
    end
    return k+1
end

rand(rng::AbstractRNG, dynsmp::DynamicSampler) = getel(dynsmp, rand(rng))
rand(dynsmp::DynamicSampler) = rand(Random.GLOBAL_RNG, dynsmp)

@inline getindex(dynsmp::DynamicSampler, i::Int) = dynsmp.v[i]

@propagate_inbounds function setindex!(dynsmp::DynamicSampler, x::Float64, i::Int)
    @extract dynsmp : v ps N levs tinds tpos
    @boundscheck 1 ≤ i ≤ N || throw(BoundsError(dynsmp, i))

    if dynsmp.trefresh ≥ max(N, 100)
        refresh!(dynsmp)
    end
    dynsmp.trefresh += 1

    @inbounds d = x - v[i]
    @inbounds v[i] = x
    dynsmp.z += d

    @inbounds @simd for j = tinds[i]:(tinds[i+1]-1)
        k = tpos[j]
        ps[k] += d
    end
    return dynsmp

    ## following is the version without using
    ## the tinds/tpos tables, for reference
    ## (may be useful to save memory)

    # k = 0
    # u = 1 << (levs-1)
    # i -= 1
    # off = 1
    # @inbounds for lev = 1:levs
    #     if i & u == 0
    #         ps[off+k] += d
    #         k *= 2
    #     else
    #         k = 2k + 1
    #     end
    #     u >>>= 1
    #     off *= 2
    # end
    # return dynsmp
end

end # module
