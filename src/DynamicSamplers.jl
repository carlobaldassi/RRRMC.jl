module DynamicSamplers

using ExtractMacro
using ..Common

import Base: getindex, setindex!, @propagate_inbounds, rand, copy!, copy

export DynamicSampler, refresh!, getel

type DynamicSampler
    v::Vec
    ps::Vec
    z::Float64
    N::Int
    levs::Int
    function DynamicSampler(N::Int)
        levs = ceil(Int, log2(N))
        N2 = 2^levs
        v2 = zeros(N2)
        ps = zeros(N2 - 1)
        return new(v2, ps, 0.0, N, levs)
    end
    function DynamicSampler(v)
        isempty(v) && throw(ArgumentError("input must be non-empty"))
        N = length(v)
        levs = ceil(Int, log2(N))
        N2 = 2^levs
        v2 = zeros(N2)
        for (i,x) in enumerate(v)
            v2[i] = x
        end
        n = findfirst(v2 .< 0)
        n == 0 || throw(ArgumentError("negative element given: v[$n] = $(v2[n])"))
        ps = zeros(N2 - 1)

        dynsmp = new(v2, ps, 0.0, N, levs)
        refresh!(dynsmp)

        return dynsmp
    end
end

# this could probably be improved a little
function refresh!(dynsmp::DynamicSampler)
    @extract dynsmp : v ps levs
    N = length(v)
    dynsmp.z = sum(v)
    L = 1
    K = N
    off = 0
    for lev = 1:levs
        @assert L < N
        i = 0
        #p = ps[lev]
        for l = off + (1:L)
            x = 0.0
            for k = 1:(K÷2)
                i += 1
                x += v[i]
            end
            i += K ÷ 2
            ps[l] = x
        end
        off += L
        L *= 2
        K ÷= 2
    end
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
    #@assert k < N
    #@assert v[k+1] > 0
    return k+1
end

rand(rng::AbstractRNG, dynsmp::DynamicSampler) = getel(dynsmp, rand(rng))
rand(dynsmp::DynamicSampler) = rand(Base.GLOBAL_RNG, dynsmp)

@inline getindex(dynsmp::DynamicSampler, i::Int) = dynsmp.v[i]

@propagate_inbounds function setindex!(dynsmp::DynamicSampler, x::Float64, i::Int)
    @extract dynsmp : v ps N levs
    @boundscheck 1 ≤ i ≤ N || throw(BoundsError(dynsmp, i))

    @inbounds d = x - v[i]
    @inbounds v[i] = x
    dynsmp.z += d
    k = 0
    u = 1 << (levs-1)
    i -= 1
    off = 1
    @inbounds for lev = 1:levs
        if i & u == 0
            ps[off+k] += d
            k *= 2
        else
            k = 2k + 1
        end
        u >>>= 1
        off *= 2
    end
    return dynsmp
end

end # module
