module RE

using ExtractMacro
using ..Interface
using ..Common

export GraphRE, GraphRepl, Renergies

import ..Interface: energy, delta_energy, neighbors, allΔE,
                    update_cache!, delta_energy_residual

import Base: start, next, done, length, eltype

const MAXDIGITS = 8

discr{ET}(::Type{ET}, x::Real) = convert(ET, round(x, MAXDIGITS))
#discr(::Type{DFloat64}, x::Real) = x
discr{ET<:Integer}(::Type{ET}, x::Integer) = convert(ET, x)

function logcoshratio(a, b)
    #return log(cosh(a) / cosh(b))
    #log((exp(a) + exp(-a)) / (exp(b) + exp(-b)))
    a = abs(a)
    b = abs(b)
    return a - b + (log1p(exp(-2a)) - log1p(exp(-2b)))
end

type GraphRE{M,γ,β} <: DiscrGraph{Float64}
    N::Int
    Nk::Int
    μ::IVec
    cache::LocalFields{Float64}

    function GraphRE(N::Integer, μ::Union{Void,IVec} = nothing)
        isa(M, Int) || throw(ArgumentError("invalid parameter M, expected Int, given: $(typeof(M))"))
        M > 2 || throw(ArgumentError("M must be greater than 2, given: $M"))
        isa(β, Float64) || throw(ArgumentError("invalid parameter β, expected Float64, given: $(typeof(β))"))
        isa(γ, Float64) || throw(ArgumentError("invalid parameter γ, expected Float64, given: $(typeof(γ))"))
        N % M == 0 || throw(ArgumentError("N must be divisible by M, given: N=$N M=$M"))
        Nk = N ÷ M

        μ::IVec = (μ ≡ nothing ? zeros(Int, Nk) : μ)
        @assert length(μ) == Nk
        cache = LocalFields{Float64}(N)
        return new(N, Nk, μ, cache)
    end
end

@doc """
    GraphRE{M,γ,β}(N::Integer) <: DiscrGraph

    TODO
""" -> GraphRE{M,γ,β}(N::Integer)

GraphRE{M,oldγ,β}(X::GraphRE{M,oldγ,β}, newγ::Float64) = GraphRE{M,newγ,β}(X.N, X.μ)

function energy{M,γ,β}(X::GraphRE{M,γ,β}, C::Config)
    @assert X.N == C.N
    @extract X : Nk μ cache
    @extract cache : lfields lfields_last
    @extract C : s

    fill!(μ, 0)
    j = 0
    for k = 1:M, i = 1:Nk
        j += 1
        @assert j == i + (k-1) * Nk
        σj = 2s[j] - 1
        μ[i] += σj
    end

    n = 0.0
    for i = 1:Nk
        n -= log(2cosh(γ * μ[i])) / β
    end

    j = 0
    for k = 1:M, i = 1:Nk
        j += 1
        @assert j == i + (k-1) * Nk
        σj = 2s[j] - 1
        μ̄ = μ[i] - σj
        kj = fk(μ̄, γ, β)
        lfields[j] = σj * kj
    end
    cache.move_last = 0
    fill!(lfields_last, 0.0)
    return n
end

function update_cache!{M,γ,β}(X::GraphRE{M,γ,β}, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract C : N s
    @extract X : Nk μ cache

    @extract cache : lfields lfields_last move_last
    if move_last == move
        @inbounds begin
            for y in neighbors(X, move)
                lfields[y], lfields_last[y] = lfields_last[y], lfields[y]
            end
            lfields[move] = -lfields[move]
            lfields_last[move] = -lfields_last[move]
        end
        σx = 2s[move] - 1
        i = mod1(move, Nk)
        μ[i] += 2σx
        return
    end

    @inbounds begin
        Ux = neighbors(X, move)
        for y in Ux
            lfields_last[y] = lfields[y]
        end
        σx = 2s[move] - 1
        i = mod1(move, Nk)
        μnew = μ[i] + 2σx
        μ[i] = μnew
        for y in Ux
            # @assert y ≠ move
            σy = 2s[y] - 1
            μ̄ = μnew - σy
            # @assert -M+1 ≤ μ̄ ≤ M-1    μ̄
            ky = fk(μ̄, γ, β)
            lfields[y] = discr(Float64, σy * ky)
        end
        lfm = lfields[move]
        lfields_last[move] = lfm
        lfields[move] = -lfm
    end
    cache.move_last = move

    # lfields_bk = copy(lfields)
    # energy(X, C)
    # lfields_bk ≠ lfields && @show move hcat(lfields,lfields_bk)
    # @assert lfields_bk == lfields

    return
end

function delta_energy{γ}(X::GraphRE{γ}, C::Config, move::Int)
    @assert X.N == C.N
    @assert 1 ≤ move ≤ C.N
    @extract X : cache
    @extract cache : lfields

    @inbounds Δ = lfields[move]
    return Δ
end

immutable CavityRange
    j0::Int
    j1::Int
    jX::Int
    st::Int
    function CavityRange(j0::Integer, j1::Integer, jX::Integer, st::Integer)
        j0 ≤ jX ≤ j1 || throw(ArgumentError("invalid CavityRange parameters, expected j0≤jX≤j1, given: j0=$j0, j1=$j1, jX=$X"))
        # TODO check step
        return new(j0, j1, jX, st)
    end
end

start(crange::CavityRange) = crange.j0 + (crange.jX == crange.j0) * crange.st
done(crange::CavityRange, j) = j > crange.j1
@inline function next(crange::CavityRange, j)
    @extract crange : j0 j1 jX st
    @assert j ≠ jX
    nj = j + st
    nj += st * (nj == jX)
    return (j, nj)
end
length(crange::CavityRange) = (crange.j1 - crange.j0) ÷ crange.st
eltype(::Type{CavityRange}) = Int

@inline function neighbors{M}(X::GraphRE{M}, j::Int)
    @extract X : Nk
    j0 = mod1(j, Nk)
    j1 = j0 + Nk * (M-1)
    return CavityRange(j0, j1, j, Nk)
end

fk(μ̄::Integer, γ::Real, β::Real) = discr(Float64, logcoshratio(γ * (μ̄ + 1), γ * (μ̄ - 1)) / β)
@generated function allΔE{M,γ,β}(::Type{GraphRE{M,γ,β}})
    K = M - 1
    iseven(K) ? Expr(:tuple, ntuple(d->fk(2*(d-1), γ, β), K÷2+1)...) :
                Expr(:tuple, ntuple(d->fk(2d-1, γ, β), (K+1)÷2)...)
end

##  # Add Transverse field to (almost) any AbstractGraph
##
##  type GraphRepl{γ,G<:AbstractGraph} <: DoubleGraph{Float64}
##      N::Int
##      M::Int
##      Nk::Int
##      X0::GraphRE{γ}
##      X1::Vector{G}
##      C1::Vector{Config}
##      β::Float64
##      Γ::Float64
##      function GraphRepl(N::Integer, M::Integer, β::Float64, Γ::Float64, g0::G, Gconstr, args...)
##          X0 = GraphRE{γ}(N, M)
##          Nk = X0.Nk
##          #J = gen_J(Nk)
##          X1 = Array{G}(M)
##          X1[1] = g0
##          for k = 2:M
##              X1[k] = Gconstr(args...)
##          end
##          C1 = [Config(Nk, init=false) for k = 1:M]
##          return new(N, M, Nk, X0, X1, C1, β, Γ)
##      end
##  end
##
##  #  """
##  #      GraphQIsingT(N::Integer, M::Integer, Γ::Float64, β::Float64) <: DoubleGraph
##  #
##  #  A `DoubleGraph` which implements a quantum Ising spin model in a transverse magnetic field,
##  #  using the Suzuki-Trotter transformation.
##  #  `N` is the number of spins, `M` the number of Suzuki-Trotter replicas, `Γ` the transverse
##  #  field, `β` the inverse temperature.
##  #  The graph is fully-connected, the interactions are random (\$J ∈ {-1,1}\$),
##  #  there are no external longitudinal fields.
##  #
##  #  """
##  function GraphRepl(Nk::Integer, M::Integer, Γ::Float64, β::Float64, Gconstr, args...)
##      @assert Γ ≥ 0
##      γ = round(2/β * log(coth(β * Γ / M)), MAXDIGITS)
##      # H0 = Nk * M / 2β * log(sinh(2β * Γ / M) / 2) # useless!!!
##      g0 = Gconstr(args...)
##      G = typeof(g0)
##      return GraphRepl{γ,G}(Nk * M, M, β, Γ, g0, Gconstr, args...)
##  end
##
##  function update_cache!(X::GraphRepl, C::Config, move::Int)
##      @extract X : X1 Nk C1
##      #@extract C : s
##      k = (move - 1) ÷ Nk + 1
##      i = mod1(move, Nk)
##
##      spinflip!(X1[k], C1[k], i)
##
##      #s1 = C1[k].s
##      #copy!(s1, 1, s, (k-1)*Nk + 1, Nk)
##      #@assert C1[k].s == s[((k-1)*Nk + 1):k*Nk]
##  end
##
##  function energy(X::GraphRepl, C::Config)
##      @assert X.N == C.N
##      @extract X : M Nk X0 X1 C1
##      @extract C : s
##
##      E = energy(X0, C)
##
##      for k = 1:M
##          s1 = C1[k].s
##          copy!(s1, 1, s, (k-1)*Nk + 1, Nk)
##          E += energy(X1[k], C1[k]) / M
##      end
##
##      return E
##  end
##
##  function Renergies(X::GraphRepl)
##      @extract X : M X1 C1
##
##      Es = zeros(M)
##
##      for k = 1:M
##          Es[k] = energy(X1[k], C1[k])
##      end
##
##      return Es
##  end
##
##  function delta_energy_residual(X::GraphRepl, C::Config, move::Int)
##      @extract X : M Nk X1 C1
##      #@extract C : s
##
##      k = (move - 1) ÷ Nk + 1
##      #s1 = C1[k].s
##      #copy!(s1, 1, s, (k-1)*Nk + 1, Nk)
##      #@assert C1[k].s == s[((k-1)*Nk + 1):k*Nk]
##
##      i = mod1(move, Nk)
##      return delta_energy(X1[k], C1[k], i) / M
##  end
##
##  function delta_energy(X::GraphRepl, C::Config, move::Int)
##      return delta_energy(X.X0, C, move) +
##             delta_energy_residual(X, C, move)
##  end

end
