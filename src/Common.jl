# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module Common

using Compat

export Vec, Vec2, IVec, unsafe_bitflip!, discretize,
       LocalFields

typealias Vec  Vector{Float64}
typealias Vec2 Vector{Vec}
typealias IVec Vector{Int}

@inline function unsafe_bitflip!(Bc::Array{UInt64}, i::Int)
    i1, i2 = Base.get_chunks_id(i)
    u = UInt64(1) << i2
    @inbounds begin
        #Bc[i1] ⊻= u
        Bc[i1] = Bc[i1] ⊻ u
    end
end
@inline unsafe_bitflip!(B::BitArray, i::Int) = unsafe_bitflip!(B.chunks, i)

type LocalFields{ET}
    lfields::Vector{ET}
    lfields_last::Vector{ET}
    move_last::Int
    function LocalFields(N::Int)
        lfields = zeros(ET, N)
        lfields_last = zeros(ET, N)
        return new(lfields, lfields_last, 0)
    end
end

function discretize{T<:Real,S<:Real}(x::T, LEV::Tuple{S,Vararg{S}})
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

function discretize{T<:Real,S<:Real}(cvec::Vector{T}, LEV::Tuple{S,Vararg{S}})
    N = length(cvec)
    fields = Array{S}(N)
    rfields = Array{T}(N)
    for (i,x) in enumerate(cvec)
        d, r = discretize(x, LEV)
        fields[i] = d
        rfields[i] = r
    end
    return fields, rfields
end

function discretize{N,T<:Real,S<:Real}(cvec::NTuple{N,T}, LEV::Tuple{S,Vararg{S}})
    fields = Array{S}(N)
    rfields = Array{T}(N)
    for (i,x) in enumerate(cvec)
        d, r = discretize(x, LEV)
        fields[i] = d
        rfields[i] = r
    end
    return tuple(fields...), tuple(rfields...)
end

end # module
