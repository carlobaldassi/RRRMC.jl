# This file is a part of RRRMC.jl. License is MIT: http://github.com/carlobaldassi/RRRMC.jl/LICENCE.md

"""
    module RRRMC

This module implements methods for reduced-rejection-rate Monte Carlo on Ising spin models.

See [`standardMC`](@ref), [`rrrMC`](@ref), [`bklMC`](@ref), [`wtmMC`](@ref) and [`extremal_opt`](@ref).
"""
module RRRMC

export standardMC, rrrMC, bklMC, wtmMC, extremal_opt

using ExtractMacro

include("DFloats.jl")
include("Common.jl")
include("Interface.jl")
include("DeltaE.jl")
include("WaitingTimes.jl")

using .Common
using .Interface
using .DeltaE
using .WaitingTimes

include("load_graphs.jl")

include("QAliases.jl")
include("REAliases.jl")
include("LEAliases.jl")
include("TLEAliases.jl")
using .QAliases
using .REAliases
using .LEAliases
using .TLEAliases

@inline accept(x) = x ≥ 0 || rand() < exp(x)
@inline function accept(c, x)
    c ≥ 1 && x ≥ 0 && return true
    a = c * exp(x)
    return a ≥ 1 || rand() < a
end

"""
    standardMC(X::AbstractGraph, β::Real, iters::Integer; keywords...)

Runs `iters` iterations of a standard Metropolis Monte Carlo algorithm on the given Ising spin model `X`, at inverse temperature `β`.
Each spin flip attempt counts as an iteration.

Returns two objects: a vector of energies and the last configuration (see [`Config`](@ref)).

Possible keyord arguments are:

* `step`: the interval of iterations with which to collect energies (for the returned result) and to call the `hook` function (see below).
  The default is `1`, which is good for debugging but otherwise generally a bad idea.
* `seed`: the random seed. The default is some arbitrary number.
* `C0`: the initial configuration. The default is `nothing`, in which case it is initialized at random. Otherwise it can be a [`Config`](@ref) object.
  Passing the result of a previous run can be useful e.g. when implementing a simulated annealing protocol, or if the system has not equilibrated yet.
* `hook`: a function to be executed after every `step` number of iterations (see above). It must take five arguments: the current iteration, the graph `X`,
  the current configuration, the number of accepted moves so far, and the current energy. Useful to collect data other than the energy, write to files ecc;
  you'd probably want to use a closure, see the example below. The return value must be a `Bool`: return `false` to interrupt the simulation, `true` otherwise.
  The default does nothing and just returns `true`.

Basic example:

```
julia> srand(76543); X = RRRMC.GraphPSpin3(3999, 5); β = 1.0;
julia> Es, C = standardMC(X, β, 100_000, step = 1_000);
```

Example of using the `hook` for collecting samples as the columns of a `BitMatrix`:

```
julia> iters = 100_000; step = 1_000; l = iters ÷ step; N = RRRMC.getN(X);
julia> Cs = BitArray(N, l); hook = (it, X, C, acc, E) -> (Cs[:,it÷step]=C.s; true);
julia> Es, C = standardMC(X, β, iters, step = step, hook = hook);
```
"""
function standardMC{ET}(X::AbstractGraph{ET}, β::Real, iters::Integer;
                        seed = 167432777111,
                        step::Integer = 1,
                        hook = (x...)->true,
                        C0::Union{Config,Void} = nothing,
                        pp = nothing,
                        quiet::Bool = false
                       )
    seed > 0 && srand(seed)
    Es = empty!(Array{ET}(min(10^8, iters ÷ step)))

    N = getN(X)
    C::Config = C0 ≡ nothing ? Config(N) : C0
    C.N == N || throw(ArgumentError("Invalid C0, wrong N, expected $N, given: $(C.N)"))
    E = energy(X, C)
    accepted = 0
    # pp ≡ nothing && (pp = randperm(N))
    # j = 0
    it = 0
    while it < iters
        it += 1
        # println("it=$it")
        # @assert isapprox(E , energy(X, C), atol=1e-8) (E, energy(X, C))
        if (it % step == 0)
            # println("it=$it")
            # @assert isapprox(E , energy(X, C), atol=1e-8) (E, energy(X, C))
            push!(Es, E)
            hook(it, X, C, accepted, E) || break
        end
        # j += 1
        # j > N && (j -= N)
        # i = pp[j]
        i = rand(1:N)
        ΔE = delta_energy(X, C, i)
        accept(-β * ΔE) || continue
        spinflip!(X, C, i)
        E += ΔE
        accepted += 1
    end
    if !quiet
        println("samples = ", length(Es))
        println("iters = ", it)
        println("accept rate = ", accepted / it)
    end
    # @show E, energy(X, C), abs(E - energy(X, C))
    return Es, C
end

### Begin RRR-related functions

function step_rrr(X, C::Config, ΔEcache)
    z = get_z(ΔEcache)
    move, ΔE = rand_move(ΔEcache)
    compute_staged!(X, C, move, ΔEcache)
    z′ = compute_reverse_probabilities!(ΔEcache)
    c = z / z′
    return move, c, ΔE
end

"""
    rrrMC(X::AbstractGraph, β::Real, iters::Integer; keywords...)

Same as [`standardMC`](@ref), but uses the reduced-rejection-rate method. Each iteration takes more time, but has a higher chance of being accepted,
so fewer iterations overall should be needed normally. Whether this trade-off is convenient depends on the parameters and the details of the model.
This function has specialized versions for [`DiscrGraph`](@ref) and [`DoubleGraph`](@ref) models.

The return values and the keyword arguments are the same as [`standardMC`](@ref), see the usage examples for that function.
"""
function rrrMC{ET}(X::SingleGraph{ET}, β::Real, iters::Integer;
                   seed = 167432777111,
                   step::Integer = 1,
                   hook = (x...)->true,
                   C0::Union{Config,Void} = nothing,
                   staged_thr::Real = NaN,
                   staged_thr_fact::Real = 5.0,
                   quiet::Bool = false
                  )

    isfinite(β) || throw(ArgumentError("β must be finite, given: $β"))
    seed > 0 && srand(seed)
    Es = empty!(Array{ET}(min(10^8, iters ÷ step)))

    if staged_thr ≡ NaN
        staged_thr = isa(X, SimpleGraph) ? 0.8 : 0.5
    end

    N = getN(X)
    C::Config = C0 ≡ nothing ? Config(N) : C0
    C.N == N || throw(ArgumentError("Invalid C0, wrong N, expected $N, given: $(C.N)"))
    E = energy(X, C)
    ΔEcache = gen_ΔEcache(X, C, β)
    #check_consistency(ΔEcache)

    λ = staged_thr_fact / N
    staged_its = 0

    it = 0
    accepted = 0
    acc_rate = 0.5
    while it < iters
        #@assert E == energy(X, C)
        #check_consistency(ΔEcache)
        it += 1
        if (it % step == 0)
            push!(Es, E)
            hook(it, X, C, accepted, E) || break
        end
        acc = false
        if acc_rate < staged_thr
            staged_its += 1
            move, c, ΔE = step_rrr(X, C, ΔEcache)
            if rand() < c
                spinflip!(X, C, move)
                apply_staged!(ΔEcache)
                E += ΔE
                accepted += 1
                acc = true
            end
        else
            move, ΔE = rand_move(ΔEcache)
            c = apply_move!(X, C, move, ΔEcache)
            if rand() < c
                E += ΔE
                accepted += 1
                acc = true
            else
                apply_move!(X, C, move, ΔEcache)
            end
        end
        acc_rate = acc_rate * (1 - λ) + acc * λ
    end
    if !quiet
        println("samples = ", length(Es))
        println("iters = ", it)
        println("accept rate = ", accepted / it)
        println("frac. staged iters = ", staged_its / it)
    end
    return Es, C
end

function rrrMC{GT,ET}(X::DoubleGraph{GT,ET}, β::Real, iters::Integer;
                      seed = 167432777111,
                      step::Integer = 1,
                      hook = (x...)->true,
                      C0::Union{Config,Void} = nothing,
                      staged_thr::Real = 0.5,
                      staged_thr_fact::Real = 5.0,
                      quiet::Bool = false
                     )
    isfinite(β) || throw(ArgumentError("β must be finite, given: $β"))
    seed > 0 && srand(seed)
    Es = empty!(Array{ET}(min(10^8, iters ÷ step)))

    N = getN(X)
    C::Config = C0 ≡ nothing ? Config(N) : C0
    C.N == N || throw(ArgumentError("Invalid C0, wrong N, expected $N, given: $(C.N)"))
    E = energy(X, C)
    @assert isfinite(E)
    X0 = inner_graph(X)
    ΔEcache = gen_ΔEcache(X0, C, β)
    #check_consistency(ΔEcache)

    λ = staged_thr_fact / N
    staged_its = 0

    it = 0
    accepted = 0
    acc_rate = 0.5
    while it < iters
        #@assert abs(E - energy(X, C)) < 1e-10 (E, energy(X, C), abs(E - energy(X,C)))
        #DeltaE.check_consistency(ΔEcache)
        it += 1
        if (it % step == 0)
            push!(Es, E)
            hook(it, X, C, accepted, E) || break
        end
        acc = false
        if acc_rate < staged_thr
            staged_its += 1
            move, c, ΔE0 = step_rrr(X0, C, ΔEcache)
            ΔE1 = delta_energy_residual(X, C, move)
            if accept(c, -β * ΔE1)
                spinflip!(X, C, move)
                apply_staged!(ΔEcache)
                E += ΔE0 + ΔE1
                accepted += 1
                acc = true
            end
        else
            move, ΔE0 = rand_move(ΔEcache)
            ΔE1 = delta_energy_residual(X, C, move)
            c = apply_move!(X, C, move, ΔEcache)
            if accept(c, -β * ΔE1)
                E += ΔE0 + ΔE1
                accepted += 1
                acc = true
            else
                apply_move!(X, C, move, ΔEcache)
            end
        end
        acc_rate = acc_rate * (1 - λ) + acc * λ
    end
    if !quiet
        println("samples = ", length(Es))
        println("iters = ", it)
        println("accept rate = ", accepted / it)
        println("frac. staged iters = ", staged_its / it)
    end
    return Es, C
end

### Begin BKL-related functions

apply_step_bkl!(X::AbstractGraph, C::Config, move::Int, ΔEcache) =
    apply_move!(X, C, move, ΔEcache, Val{false})

apply_step_bkl!(X::DiscrGraph, C::Config, move::Int, ΔEcache) =
    apply_move!(X, C, move, ΔEcache)

"""
    bklMC(X::AbstractGraph, β::Real, iters::Integer; keywords...)

Same as [`standardMC`](@ref), but uses the rejection-free method by Bortz, Kalos and Lebowitz. Each step takes more time, but rejected moves
are almost free, since they are skipped entirely. This function has a specialized version for [`DiscrGraph`](@ref) models.

The return values and the keyword arguments are the same as [`standardMC`](@ref), see the usage examples for that function.

Note that the number of iterations includes the rejected moves. This makes the results directly comparable with those of `standardMC`. It also
means that increasing `β` at fixed `iters` will result in fewer steps being actually computed.
"""
function bklMC{ET}(X::AbstractGraph{ET}, β::Real, iters::Integer;
                   seed = 167432777111,
                   step::Integer = 1,
                   hook = (x...)->true,
                   C0::Union{Config,Void} = nothing,
                   quiet::Bool = false
                  )
    seed > 0 && srand(seed)
    Es = empty!(Array{ET}(min(10^8, iters ÷ step)))

    N = getN(X)
    C::Config = C0 ≡ nothing ? Config(N) : C0
    C.N == N || throw(ArgumentError("Invalid C0, wrong N, expected $N, given: $(C.N)"))
    E = energy(X, C)
    ΔEcache = gen_ΔEcache(X, C, β, false)
    #check_consistency(ΔEcache)

    it = 0
    accepted = 0
    nextstep = step
    stop = false
    while it < iters
        #@assert E == energy(X, C)
        #check_consistency(ΔEcache)

        skip = rand_skip(ΔEcache)
        move, ΔE = rand_move(ΔEcache)

        while it + skip + 1 ≥ nextstep
            push!(Es, E)
            hook(nextstep, X, C, accepted, E) || @goto out
            nextstep += step
            nextstep > iters && @goto out
        end

        apply_step_bkl!(X, C, move, ΔEcache)
        it += skip + 1
        E += ΔE
        accepted += 1
    end
    @label out
    if !quiet
        println("samples = ", length(Es))
        println("iters = ", it)
        println("accept rate = ", accepted / iters)
        println("true it = ", accepted)
    end
    return Es, C
end

"""
    wtmMC(X::AbstractGraph, β::Real, samples::Integer; keywords...)

Same as [`standardMC`](@ref), but uses the rejection-free waiting-time method by Dall and Sibani. It is similar to [`bklMC`](@ref).

The return values and the keyword arguments are *almost* the same as [`standardMC`](@ref), see the usage examples for that function.
However, the waiting time method uses an internal "global time" variable, which takes the place of the "iterations" counter of [`standardMC`](@ref)
and of [`bklMC`](@ref). Thus, this function has two differences with respect to the other samplers in the module:

* the function takes a `samples` integer argument, with the maximum number of samples which will be collected.
* the `step` keyword argument is of type `Float64` instead of integer (default value = `1.0`; you'll probably want to change this).
  The `step` is measured in terms of the global time variable and is scaled with the size of the problem `N`.

The total number of samples actually collected can still be less than `samples` if the `hook` function from the keyword arguments returns `false` earlier.
"""
function wtmMC{ET}(X::AbstractGraph{ET}, β::Real, samples::Integer;
                   seed = 167432777111,
                   step::Float64 = 1.0,
                   hook = (x...)->true,
                   C0::Union{Config,Void} = nothing,
                   quiet::Bool = false
                  )
    seed > 0 && srand(seed)
    Es = empty!(Array{ET}(min(10^8, samples)))

    N = getN(X)
    C::Config = C0 ≡ nothing ? Config(N) : C0
    C.N == N || throw(ArgumentError("Invalid C0, wrong N, expected $N, given: $(C.N)"))
    E = energy(X, C)
    theap = THeap(X, C, β)

    step /= N

    num_moves = 0
    tmax = step * samples
    t = 0.0
    nextstep = step
    while t < tmax
        #@assert abs(E - energy(X, C)) < 1e-10 E-energy(X,C)

        t′, move = pick_next(theap)
        while t′ ≥ nextstep
            push!(Es, E)
            hook(nextstep, X, C, num_moves, E) || @goto out
            nextstep += step
            nextstep > tmax + 1e-10 && @goto out
        end
        t = t′
        ΔE = update_heap!(theap, X, C, β, move, t)
        E += ΔE

        num_moves += 1
    end
    @label out
    if !quiet
        println("samples = ", length(Es))
        println("num_moves = ", num_moves)
        println("global time = ", t)
        println("ratio = ", t / num_moves)
    end
    return Es, C
end



# τ-Extremal Optimization

"""
    extremal_opt(X::AbstractGraph, τ::Real, iters::Integer; keywords...)

Extremal optimization algorithm: it seeks the lowest energy state of a given Ising spin model `X` by performing
a random walk biased towards low-energy states, but with a heavy tail which allows it to avoid getting trapped.

The interface is very similar to that of [`standardMC`](@ref) and the other RRRMC Monte Carlo functions, but it
takes a parameter `τ` instead of `β`, the return values are different, and the `hook` keyword argument has a
different signature, see below.

The parameter `τ` controls the shape of the tail and should be larger than 1 (a reasonable value could be 1.3);
`iters` is the total number of spin flips performed.

This function has a specialized version for [`DiscrGraph`](@ref) models; otherwise, it works but it is not
implemented very efficiently at the moment (each spin flip takes O(N) time even on diluted graphs).

Returns 4 objects: the final configuration, the minimum energy found, the configuration of minimum energy, and the
iteration at which such configuration was found (see [`Config`](@ref)).

Possible keyord arguments are:
* `step`: the interval of iterations with which to call the `hook` function (see below).
  The default is `1`, which is good for debugging but otherwise generally a bad idea, probably.
* `seed`: the random seed. The default is some arbitrary number.
* `C0`: the initial configuration. The default is `nothing`, in which case it is initialized at random. Otherwise it can be a [`Config`](@ref) object.
* `hook`: a function to be executed after every `step` number of iterations (see above). It must take five arguments: the current iteration, the graph `X`,
  the current configuration, the current energy, and the minimum energy found so far. Useful to collect data other than the energy, write to files ecc;
  you'd probably want to use a closure, see the example below. The return value must be a `Bool`: return `false` to interrupt the simulation, `true` otherwise.
  The default does nothing and just returns `true`. Note that the signature is similar to that of the corresponding `hook` argument in [`standardMC`](@ref),
  but different.

Basic example:

```
julia> srand(76543); X = RRRMC.GraphPSpin3(3999, 5); τ = 1.3;
julia> C, Emin, Cmin, itmin = extremal_opt(X, τ, 100_000, step = 1_000);
```

Example of using the `hook` for collecting samples as the columns of a `BitMatrix`:

```
julia> iters = 100_000; step = 1_000; l = iters ÷ step; N = RRRMC.getN(X);
julia> Cs = BitArray(N, l); hook = (it, X, C, E, Emin) -> (Cs[:,it÷step]=C.s; true);
julia> C, Emin, Cmin, itmin = extremal_opt(X, τ, iters, step = step, hook = hook);
```

"""
function extremal_opt{ET}(X::AbstractGraph{ET}, τ::Real, iters::Integer;
                         seed = 167432777111,
                         step::Integer = 1,
                         hook = (x...)->true,
                         C0::Union{Config,Void} = nothing,
                         quiet::Bool = false
                        )
    seed > 0 && srand(seed)

    N = getN(X)
    C::Config = C0 ≡ nothing ? Config(N) : C0
    C.N == N || throw(ArgumentError("Invalid C0, wrong N, expected $N, given: $(C.N)"))
    E = energy(X, C)
    Emin = E
    Cmin = copy(C)
    itmin = 0
    ΔEcache = gen_EOcache(X, C, τ)
    # check_consistency(ΔEcache)

    it = 0
    stop = false
    while it < iters
        it += 1
        # @assert E == energy(X, C)
        # check_consistency(ΔEcache)

        if (it % step == 0)
            hook(it, X, C, E, Emin) || break
            # E = energy(X, C)
        end

        move, ΔE = rand_move(ΔEcache)
        apply_move!(X, C, move, ΔEcache)
        E += ΔE

        if E < Emin
            Emin = E
            copy!(Cmin, C)
            itmin = it
        end
    end
    if !quiet
        # println("samples = ", length(Es))
        println("iters = ", it)
        println("min [it = $itmin] = $Emin")
    end
    return C, Emin, Cmin, itmin
end

### auxiliary/miscellanaous/ugly stuff. Hic sunt leones.

ba2int(b::BitArray) = Int(b.chunks[1])
int2ba(i::Int, N) = (@assert N≤64; b = BitVector(N); b.chunks[1]=i; b)

function truep(X::DiscrGraph, β::Real)
    N = X.N
    @assert N ≤ 64
    S = 1 << N
    C = Config(N)
    p = zeros(S)
    Z = 0.0
    for i = 1:S
        C.s.chunks[1] = i - 1
        E = energy(X, C)
        p[i] = exp(-β * E)
        Z += p[i]
    end
    scale!(p, 1/Z)
    return p
end

function samplep(ws::Vector{BitVector})
    N = length(ws[1])
    S = 1 << N
    C = MCNew.Config(N)
    p = zeros(S)
    for w in ws
        copy!(C.s, w)
        i = ba2int(w) + 1
        p[i] += 1
    end
    scale!(p, 1/length(ws))
    return p
end

function tm{T<:Real}(Es::Vector{T}, step = 1, skip0 = 0.1, skip1 = 0.05)
    N = length(Es)
    i0 = floor(Int, N * skip0)
    N0 = N - i0
    n = N0 ÷ step
    m = zeros(n)
    m[1] = mean(Es[i0 + (1:step)])
    for j = 2:n
        m[j] = (mean(Es[i0 + (j-1)*step + (1:step)]) + m[j-1] * (j-1)) / j;
    end
    return m[1 + floor(Int, skip1 * n):end]
end

function ravg{T<:Real}(Es::Vector{T}, step = 1_000, skip0 = 0.0)
    N = length(Es)
    i0 = floor(Int, N * skip0)
    N0 = N - i0
    n = N0 ÷ step
    m = zeros(n)
    m[1] = mean(Es[i0 + (1:step)])
    for j = 2:n
        m[j] = mean(Es[i0 + (j-1)*step + (1:step)])
    end
    return m
end

function second_eigenvalue(Q::AbstractMatrix{Float64})
    ev = eigvals(Q)
    @assert all(x->(abs(imag(x)) ≤ 1e-10), ev)
    evr = Float64[real(x) for x in ev]

    τ = -1 / log(sort!(evr)[end-1])
end

function second_eigenvalue_standard(X::DiscrGraph, β::Float64)
    N = getN(X)
    @assert N ≤ 64
    S = 1 << N
    Q = zeros(S, S)
    C = MCNew.Config(N)
    p = zeros(S)
    for ks = 1:S
        C.s.chunks[1] = ks - 1
        pT = 0.0
        for j = 1:N
            kd = Int(((ks-1) $ (UInt64(1) << (j-1))) + 1)

            #=
            C.s[j] $= 1
            @assert C.s.chunks[1] == kd
            C.s[j] $= 1
            kd += 1
            =#

            ΔE = delta_energy(X, C, j)
            p = min(1.0, exp(-β * ΔE)) / N
            pT += p
            Q[kd, ks] = p
        end
        Q[ks, ks] = 1 - pT
    end

    τ = second_eigenvalue(Q)

    return Q, τ
end

function second_eigenvalue_bkl(Q::Matrix{Float64})
    pr = diag(Q)
    rfQ = (Q - diagm(pr)) ./ (1 - pr')
    τ = second_eigenvalue(rfQ)

    return rfQ, τ
end

function second_eigenvalue_bkl(X::DiscrGraph, β::Float64)
    Q, _ = second_eigenvalue_standard(X, β)
    return second_eigenvalue_bkl(Q)
end

function second_eigenvalue_rrr(X::DiscrGraph, β::Float64)
    N = getN(X)
    @assert N ≤ 64
    S = 1 << N

    Q = zeros(S, S)
    C = Config(N)

    for ks = 1:S
        C.s.chunks[1] = ks - 1

        energy(X, C)
        ΔEcache = gen_ΔEcache(X, C, β)
        #DeltaE.check_consistency(ΔEcache)

        pchg = 0.0
        for j = 1:N
            @extract ΔEcache : k=pos[j] z
            compute_staged!(X, C, j, ΔEcache)
            z′ = compute_reverse_probabilities!(ΔEcache)
            w = get_class_f(ΔEcache, k)

            # pp = min(1.0, z / z′) * w / z # this is error-prone when z is small
            pp = w / max(z, z′)
            @assert isfinite(pp)
            pchg += pp

            kd = Int(((ks-1) $ (UInt64(1) << (j-1))) + 1)
            Q[kd, ks] = pp
        end
        @assert -1e-15 ≤ 1 - pchg ≤ 1 + 1e-15
        Q[ks, ks] = clamp(1 - pchg, 0.0, 1.0)
    end

    τ = second_eigenvalue(Q)

    return Q, τ
end

function second_eigenvalue_stats(;seed::Integer = 86823, graph = GraphRRG, args = (10, 5), βr = 0.5:0.5:3.0, n = 10)
    τs = Matrix{Float64}[]
    rrs = Matrix{Float64}[]
    for j = 1:n
        srand(seed + j)
        println("seed = $(seed + j)")
        X = graph(args...)
        τ = Array{Float64}(length(βr), 3)
        rr = Array{Float64}(length(βr), 3)
        for (l,β) in enumerate(βr)
            println("  β=$β")
            p = truep(X, β);
            Q, τ[l,1] = second_eigenvalue_standard(X, β)
            @assert maximum(abs(p - Q * p)) < 1e-13
            rr[l,1] = sum(diag(Q) .* p)
            rfQ, τ[l,2] = second_eigenvalue_bkl(Q)
            pr = diag(Q)
            @assert maximum(abs(p .* (1 - pr) - rfQ * (p .* (1 - pr)))) < 1e-13
            rr[l,2] = 0.0
            Q, τ[l,3] = second_eigenvalue_rrr(X, β)
            @assert maximum(abs(p - Q * p)) < 1e-13
            rr[l,3] = sum(diag(Q) .* p)
            println("    τ:  ", τ[l,:])
            println("    τr: ", τ[l,1] ./ τ[l,:])
            println("    rr: ", rr[l,:])
        end
        push!(τs, τ)
        push!(rrs, rr)
    end
    println("---------\n")
    println("τ stats1")
    τm = mean(τs)
    for (l,β) in enumerate(βr)
        println("  β=$β")
        println("    ", τm[l,:])
        println("    ", τm[l,1] ./ τm[l,:])
    end
    println("---------\n")
    println("τ stats2")
    τr = mean(Matrix{Float64}[τ[:,1] ./ τ for τ in τs])
    for (l,β) in enumerate(βr)
        println("  β=$β")
        println("    ", τr[l,:])
    end
    println("---------\n")
    println("rr stats")
    rrm = mean(rrs)
    for (l,β) in enumerate(βr)
        println("  β=$β")
        println("    ", rrm[l,:])
    end
    println("---------\n")
    return τs, rrs
end

function runtest(;N = 1_000, D = 4, β = 1.0, seedlist = 101:110, iters = 10_000_000, stepfact = 10, step = 100)
    srand(111)
    X = GraphRRG(N, D)

    nsamples = iters ÷ step

    nseeds = length(seedlist)
    EsBKL = zeros(nsamples, nseeds)
    Es5 = zeros(nsamples, nseeds)
    elbkl = 0.0
    elsmrt = 0.0
    for (l,seed) in enumerate(seedlist)
        info("SEED = $seed")
        #info("  bkl")
        #el = @elapsed EsBKL[:,l] = bkl(X, β, iters * stepfact, seed=seed, step=step*stepfact)
        #info("    time: $el")
        #elbkl += el
        info("  rrr")
        el = @elapsed Es5[:,l] = rrrMC(X, β, iters, seed=seed, step=step)
        info("    time: $el")
        elsmrt += el
        println()
    end
    info("average elapsed time bkl:  ", elbkl / nseeds)
    info("average elapsed time smrt: ", elsmrt / nseeds)
    info("ratio: ", elsmrt / elbkl)
    return EsBKL, Es5
end

end
