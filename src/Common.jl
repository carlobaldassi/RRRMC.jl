# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module Common

export Vec, Vec2, IVec, IVec2, unsafe_bitflip!, flipbits!, discretize,
       LocalFields, AllButOne

import Base: iterate

const Vec  = Vector{Float64}
const Vec2 = Vector{Vec}
const IVec = Vector{Int}
const IVec2 = Vector{IVec}

@inline function unsafe_bitflip!(Bc::Array{UInt64}, i::Int)
    i1, i2 = Base.get_chunks_id(i)
    u = UInt64(1) << i2
    @inbounds begin
        #Bc[i1] ⊻= u
        Bc[i1] = Bc[i1] ⊻ u
    end
end
@inline unsafe_bitflip!(B::BitArray, i::Int) = unsafe_bitflip!(B.chunks, i)

flipbits!(B::BitArray) = B .= .!(B)

mutable struct LocalFields{ET}
    lfields::Vector{ET}
    lfields_last::Vector{ET}
    move_last::Int
    function LocalFields{ET}(N::Int) where {ET}
        lfields = zeros(ET, N)
        lfields_last = zeros(ET, N)
        return new{ET}(lfields, lfields_last, 0)
    end
end

function discretize(x::T, LEV::Tuple{S,Vararg{S}}) where {T<:Real,S<:Real}
    d = LEV[1]
    r = x - d
    for l = 2:length(LEV)
        d1 = LEV[l]
        r1 = x - d1
        if abs(r1) < abs(r)
            d, r = d1, r1
        end
    end
    return d, r
end

function discretize(cvec::Vector{T}, LEV::Tuple{S,Vararg{S}}) where {T<:Real,S<:Real}
    N = length(cvec)
    fields = Array{S}(undef, N)
    rfields = Array{T}(undef, N)
    for (i,x) in enumerate(cvec)
        d, r = discretize(x, LEV)
        fields[i] = d
        rfields[i] = r
    end
    return fields, rfields
end

function discretize(cvec::NTuple{N,T}, LEV::Tuple{S,Vararg{S}}) where {N,T<:Real,S<:Real}
    fields = Array{S}(undef, N)
    rfields = Array{T}(undef, N)
    for (i,x) in enumerate(cvec)
        d, r = discretize(x, LEV)
        fields[i] = d
        rfields[i] = r
    end
    return tuple(fields...), tuple(rfields...)
end


# iterate over 1:N except i, i.e.
#   1, 2, 3, ..., i-1, i+1, ..., N-1, N
# useful for fully-connected models
mutable struct AllButOne
    N::Int
    i::Int
    function AllButOne(N::Integer, i::Integer)
        N ≥ 1 || throw(ArgumentError("N must me larger than 1, given: $N"))
        1 ≤ i ≤ N || throw(ArgumentError("expected 1 ≤ i ≤ N, given i=$i N=$N"))
        return new(N, i)
    end
end

function iterate(n::AllButOne, j = 1 + (n.i == 1))
    j == n.N + 1 && return nothing
    return j, j + 1 + (j == n.i - 1)
end
Base.length(n::AllButOne) = n.N - 1

end # module
