# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

# A Real type which is actually an integer in disguise.
# Similar to FixedPoint, but more efficient for
# the purposes of this code.

module DFloats

export DFloat64, MAXDIGITS

const MAXDIGITS = 5

primitive type DFloat64 <: Real 64 end
const dfact = 10^MAXDIGITS

import Base: convert, ==, <, <=, *, /, +, -, round, typemin, show, promote_rule, decompose,
             zero, signbit, abs

i2d(x::Int64) = reinterpret(DFloat64, x)
d2i(x::DFloat64) = reinterpret(Int64, x)
convert(::Type{DFloat64}, x::DFloat64) = x
convert(::Type{DFloat64}, x::Real) = i2d(round(Int64, x * dfact))
convert(::Type{Bool}, x::DFloat64) = d2i(x) > 0
convert(::Type{Integer}, x::DFloat64) = Int64(d2i(x) / dfact)
convert(::Type{T}, x::DFloat64) where {T<:Real} = T(d2i(x) / dfact)
convert(::Type{DFloat64}, x::Integer) = i2d(Int64(x * dfact))
DFloat64(x::Number) = convert(DFloat64, x)

-(x::DFloat64) = i2d(-d2i(x))
-(x::DFloat64, y::DFloat64) = i2d(d2i(x) - d2i(y))
+(x::DFloat64, y::DFloat64) = i2d(d2i(x) + d2i(y))
*(a::Bool, b::DFloat64) = i2d(a * d2i(b))
*(a::DFloat64, b::Bool) = b * a
*(a::Integer, b::DFloat64) = i2d(Int64(a) * d2i(b))
*(a::DFloat64, b::Integer) = b * a
/(a::DFloat64, b::Integer) = i2d(d2i(a) ÷ Int64(b))
==(a::DFloat64, b::DFloat64) = d2i(a) == d2i(b)
<(a::DFloat64, b::DFloat64) = d2i(a) < d2i(b)
<=(a::DFloat64, b::DFloat64) = d2i(a) <= d2i(b)

zero(::Type{DFloat64}) = i2d(Int64(0))

promote_rule(::Type{Float64}, ::Type{DFloat64}) = Float64
promote_rule(::Type{DFloat64}, ::Type{Float64}) = Float64
promote_rule(::Type{Int}, ::Type{DFloat64}) = DFloat64
promote_rule(::Type{DFloat64}, ::Type{Int}) = DFloat64

round(x::DFloat64, digits::Integer) = (@assert digits == MAXDIGITS; x)
typemin(::Type{DFloat64}) = i2d(typemin(Int64))

signbit(x::DFloat64) = signbit(d2i(x))
abs(x::DFloat64) = i2d(abs(d2i(x)))

show(io::IO, x::DFloat64) = show(io, Float64(x))

# used for hashing
function decompose(x::DFloat64)
    i = d2i(x)
    m = gcd(i, dfact)
    return i ÷ m, 1, dfact ÷ m
end

end # module
