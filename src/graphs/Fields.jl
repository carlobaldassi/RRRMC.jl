module Fields

using ExtractMacro
using ..Interface
using ..Common

export GraphFields, GraphFieldsCont

import ..Interface: energy, delta_energy, neighbors, allΔE, delta_energy_residual

type GraphFields{ET,LEV} <: DiscrGraph{ET}
    N::Int
    fields::Vector{ET}
    function GraphFields(fields::Vector{ET})
        isa(LEV, Tuple{Real,Vararg{Real}}) || throw(ArgumentError("invalid level spec, expected a Tuple of Reals, given: $LEV"))
        all(f->f ∈ LEV, fields) || throw(ArgumentError("invalid field value, expected $LEV, given: $(findfirst(f->f ∉ LEV, fields))"))
        length(unique(LEV)) == length(LEV) || throw(ArgumentError("repeated levels in LEV: $LEV"))
        N = length(fields)
        return new(N, fields)
    end
end

"""
    GraphFields(N::Integer, LEV::Tuple = (1,)) <: DiscrGraph

A simple `DiscrGraph` type with `N` non-interacting variables, each of which is
subject to a local field. The fields are extracted at random from `LEV`, which
must be a `Tuple` of `Real`s.

Mostly useful for testing/debugging purposes.
"""
function GraphFields{ET<:Real}(N::Integer, LEV::Tuple{ET,Vararg{ET}})
    return GraphFields{ET,LEV}(rand(collect(LEV), N))
end

GraphFields(N::Integer) = GraphFields{Int,(1,)}(ones(Int, N))

function energy(X::GraphFields{Int,(1,)}, C::Config)
    @extract C : N s
    return N - 2 * sum(s)
end

function energy{ET}(X::GraphFields{ET}, C::Config)
    @assert X.N == C.N
    @extract X : fields
    @extract C : N s

    local E::ET = zero(ET)
    for i = 1:N
        E -= fields[i] * (2s[i] - 1)
    end

    return E
end

function delta_energy(X::GraphFields{Int,(1,)}, C::Config, move::Int)
    @assert 1 ≤ move ≤ C.N
    @extract C : N s

    return 2 * (2 * s[move] - 1)
end

function delta_energy{ET}(X::GraphFields{ET}, C::Config, move::Int)
    @assert 1 ≤ move ≤ C.N
    @assert X.N == C.N
    @extract X : fields
    @extract C : N s

    return convert(ET, 2 * fields[move] * (2s[move] - 1))
end


neighbors(X::GraphFields, i::Int) = return ()
#allΔE(::Type{GraphFields{Int,(1,)}}) = (2,)

@generated function allΔE{ET,LEV}(::Type{GraphFields{ET,LEV}})
    absLEV = sort!(unique(map(x->convert(ET, 2*abs(x)), LEV)))
    return Expr(:tuple, absLEV...)
end


type GraphFieldsCont{ET,LEV} <: DoubleGraph{Float64}
    N::Int
    X0::GraphFields{ET,LEV}
    rfields::Vec
    function GraphFieldsCont(N::Integer)
        cfields = randn(N)
        fields, rfields = discretize(cfields, LEV)
        X0 = GraphFields{ET,LEV}(fields)
        return new(N, X0, rfields)
    end
end

"""
    GraphFieldsCont(N::Integer, LEV::Tuple) <: DoubleGraph

A simple `DoubleGraph` type with `N` non-interacting variables, each of which is
subject to a local field. The fields are extracted independently from a normal
distribution with unit variance, and then are discretized using the values in `LEV`, which
must be a `Tuple` of `Real`s.

Mostly useful for testing/debugging purposes.
"""
GraphFieldsCont{ET<:Real}(N::Integer, LEV::Tuple{ET,Vararg{ET}}) = GraphFieldsCont{ET,LEV}(N)

function energy(X::GraphFieldsCont, C::Config)
    @assert X.N == C.N
    @extract X : X0 rfields
    @extract C : N s

    E0 = energy(X0, C)

    E1 = 0.0
    for i = 1:N
        E1 -= rfields[i] * (2s[i] - 1)
    end

    return convert(Float64, E0 + E1)
end

function delta_energy_residual(X::GraphFieldsCont, C::Config, move::Int)
    @assert 1 ≤ move ≤ C.N
    @assert X.N == C.N
    @extract X : rfields
    @extract C : N s

    return 2 * rfields[move] * (2s[move] - 1)
end

function delta_energy(X::GraphFieldsCont, C::Config, move::Int)
    ΔE0 = delta_energy(X.X0, C, move)
    ΔE1 = delta_energy_residual(X, C, move)
    return convert(Float64, ΔE0 + ΔE1)
end

end # module
