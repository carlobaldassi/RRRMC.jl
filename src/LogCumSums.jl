module LogCumSums

using ExtractMacro
using ..Common

import Base: getindex, setindex!, @propagate_inbounds, rand, copy!, copy

export LogCumSum, refresh!, getel

type LogCumSum
    v::Vec
    ps::Vec2
    z::Float64
    N::Int
    function LogCumSum(N::Int)
        levs = ceil(Int, log2(N))
        N2 = 2^levs
        v2 = zeros(N2)
        ps = [zeros(2^lev) for lev = 1:levs]
        return new(v2, ps, 0.0, N)
    end
    function LogCumSum(v)
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
        ps = [zeros(2^lev) for lev = 1:levs]

        lcs = new(v2, ps, 0.0, N)
        refresh!(lcs)

        return lcs
    end
end

# this could probably be improved a little
function refresh!(lcs::LogCumSum)
    @extract lcs : v ps
    N = length(v)
    levs = length(ps)
    lcs.z = sum(v)
    L = 1
    K = N
    for lev = 1:levs
        @assert L < N
        i = 0
        p = ps[lev]
        for l = 1:L
            x = 0.0
            for k = 1:(K÷2)
                i += 1
                x += v[i]
            end
            i += K ÷ 2
            p[l] = x
        end
        L *= 2
        K ÷= 2
    end
    return lcs
end

function copy!(dest::LogCumSum, src::LogCumSum)
    @extract src : v ps z N
    dest.N == N || throw(ArgumentError("dest N=$(dest.N), src N=$N"))

    copy!(dest.v, v)
    for (dp1, sp1) in zip(dest.ps, ps)
        copy!(dp1, sp1)
    end
    dest.z = src.z
    return dest
end

copy(lcs::LogCumSum) = copy!(LogCumSum(lcs.N), lcs)

# for checking purposes
function getel_naive(v::Vec, z::Float64, x::Float64)
    @assert 0 ≤ x < 1
    #@extract lcs : v z N
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

getel_naive(lcs::LogCumSum, x::Float64) = getel_naive(lcs.v, lcs.z, x)

function getel(lcs::LogCumSum, x::Float64)
    @assert 0 ≤ x < 1
    @extract lcs : v ps z N
    x *= z

    k = 0
    @inbounds for p1 in ps
        p = p1[k+1]
        k *= 2
        if x > p
            x -= p
            k += 1
        end
    end
    @assert k < N
    @assert v[k+1] > 0
    return k+1
end

rand(rng::AbstractRNG, lcs::LogCumSum) = getel(lcs, rand(rng))
rand(lcs::LogCumSum) = rand(Base.GLOBAL_RNG, lcs)

@inline getindex(lcs::LogCumSum, i::Int) = lcs.v[i]

@propagate_inbounds function setindex!(lcs::LogCumSum, x::Float64, i::Int)
    @extract lcs : v ps N
    @boundscheck 1 ≤ i ≤ N || throw(BoundsError(lcs, i))

    @inbounds d = x - v[i]
    @inbounds v[i] = x
    lcs.z += d
    k = 0
    u = 1 << (length(ps)-1)
    i -= 1
    @inbounds for p1 in ps
        if i & u == 0
            p1[k+1] += d
            k *= 2
        else
            k = 2k + 1
        end
        u >>>= 1
    end
    return lcs
end

end # module
