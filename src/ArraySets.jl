# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

# This implements an unsorted set which allows for
# fast insertion-deletion and randomly choosing one
# element.

module ArraySets

using ExtractMacro
using ..Common

export ArraySet, check_consistency

import Base: push!, delete!, rand, length, empty!,
             start, next, done, in, show

type ArraySet
    N::Int
    v::IVec
    pos::IVec
    t::Int
    ArraySet(N::Integer) = new(N, zeros(Int, N), zeros(Int, N), 0)
end

function check_consistency(aset::ArraySet)
    @extract aset : N v pos t
    @assert 0 ≤ t ≤ N
    c = 0
    for i = 1:N
        pos[i] == 0 && continue
        c += 1
        @assert 1 ≤ pos[i] ≤ t (i, pos[i], t)
        @assert v[pos[i]] == i (v[1:t], pos, pos[i], i)
    end
    @assert c == t
    for i = 1:t
        @assert v[i] ≠ 0
        @assert pos[v[i]] == i
    end
end

length(aset::ArraySet) = aset.t

in(aset::ArraySet, i) = pos[i] ≠ 0

function show(io::IO, aset::ArraySet)
    @extract aset : v t
    print(io, "ArraySet(")
    for i = 1:(t-1)
        print(io, v[i], ',')
    end
    t > 0 && print(io, v[t])
    print(io, ')')
end

function push!(aset::ArraySet, i::Integer)
    @extract aset : v pos t
    #@assert pos[i] == 0
    @inbounds begin
    t += 1
    v[t] = i
    pos[i] = t
    aset.t = t
    end
end
function delete!(aset::ArraySet, i::Integer)
    @extract aset : v pos t
    @inbounds begin
    p = pos[i]
    #@assert 0 < p ≤ t
    v[p] = v[t]
    pos[v[p]] = p
    pos[i] = 0
    aset.t -= 1
    end
end

rand(aset::ArraySet) = rand(Base.Random.globalRNG(), aset)
function rand(r::AbstractRNG, aset::ArraySet)
    @extract aset : v t
    @inbounds i = v[rand(r, 1:t)]
    return i
end

function empty!(aset::ArraySet)
    @extract aset : pos t
    fill!(pos, 0)
    aset.t = 0
    return aset
end

start(aset::ArraySet) = 1
done(aset::ArraySet, i) = i == aset.t + 1
@inline function next(aset::ArraySet, i)
    @extract aset : v
    @boundscheck 1 ≤ i ≤ aset.t
    @inbounds x = v[i]
    return (x, i + 1)
end

end # module
