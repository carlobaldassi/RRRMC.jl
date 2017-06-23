# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

module Common

using Compat

export Vec, Vec2, IVec, IVec2, unsafe_bitflip!, discretize,
       LocalFields, @inner, AllButOne

wrapin(head::Symbol, fn, T::Symbol) = Expr(:curly, fn, T)
function wrapin(head::Symbol, fn, T::Expr)
    @assert Base.Meta.isexpr(T, [:tuple, :cell1d])
    Expr(head, fn, T.args...)
end

# horrible macro to keep compatibility with both julia 0.5 and 0.6,
# while avoiding some even more horrible syntax
macro inner(T, ex)
    VERSION < v"0.6.0-dev.2643" && return esc(ex)
    @assert Base.Meta.isexpr(ex, [:(=), :function])
    @assert length(ex.args) == 2
    @assert isa(ex.args[1], Expr) && ex.args[1].head == :call
    @assert isa(ex.args[1].args[1], Symbol)
    fn = wrapin(:curly, ex.args[1].args[1], T)
    fargs = ex.args[1].args[2:end]
    body = ex.args[2]

    return esc(Expr(ex.head, wrapin(:where, Expr(:call, fn, fargs...), T), body))
end

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

type LocalFields{ET}
    lfields::Vector{ET}
    lfields_last::Vector{ET}
    move_last::Int
    @inner {ET} function LocalFields(N::Int)
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


# iterate over 1:N except i, i.e.
#   1, 2, 3, ..., i-1, i+1, ..., N-1, N
# useful for fully-connected models
type AllButOne
    N::Int
    i::Int
    function AllButOne(N::Integer, i::Integer)
        N ≥ 1 || throw(ArgumentError("N must me larger than 1, given: $N"))
        1 ≤ i ≤ N || throw(ArgumentError("expected 1 ≤ i ≤ N, given i=$i N=$N"))
        return new(N, i)
    end
end

Base.start(n::AllButOne) = 1 + (n.i == 1)
Base.done(n::AllButOne, j::Int) = j == n.N+1
Base.next(n::AllButOne, j::Int) = j, j + 1 + (j == n.i - 1)
Base.length(n::AllButOne) = n.N - 1

end # module
