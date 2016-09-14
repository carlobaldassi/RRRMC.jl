"""
    module RRRMC

This module implements methods for reduced-rejection-rate Monte Carlo on Ising spin models.

See [`standardMC`](@ref), [`rrrMC`](@ref) and [`bklMC`](@ref).
"""
module RRRMC

export standardMC, rrrMC, bklMC, wtmMC

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
using .QAliases
using .REAliases

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
  The default is a no-op and just returns `true`.

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
function standardMC{ET}(X::AbstractGraph{ET}, β::Real, iters::Integer; seed = 167432777111, step::Integer = 1, hook = (x...)->true, C0::Union{Config,Void} = nothing, pp = nothing)
    seed > 0 && srand(seed)
    Es = empty!(Array(ET, min(10^8, iters ÷ step)))

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
        #println("it=$it")
        #@assert abs(E - energy(X, C)) < 1e-8 (E, energy(X, C))
        if (it % step == 0)
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
    println("samples = ", length(Es))
    println("iters = ", it)
    println("accept rate = ", accepted / it)
    return Es, C
end

### Begin RRR-related functions

function step_rrr(X::DiscrGraph, C::Config, ΔEcache::DeltaECache)
    @extract ΔEcache : z
    move, ΔE = rand_move(ΔEcache)
    compute_staged!(X, C, move, ΔEcache)
    z′ = compute_reverse_probabilities!(ΔEcache)
    c = z / z′
    return move, c, ΔE
end

"""
    rrrMC(X::AbstractGraph, β::Real, iters::Integer; keywords...)

Same as [`standardMC`](@ref), but uses the reduced-rejection-rate method. Each iteration takes moretime, but has a higher chance of being accepted,
so fewer iterations overall should be needed normally. Whether this trade-off is convenient depends on the parameters and the details of the model.

The return values and the keyword arguments are the same as [`standardMC`](@ref), see the usage examples for that function. Note however
that this function can only be used with [`DiscrGraph`](@ref) or [`DoubleGraph`](@ref) models.
"""
function rrrMC{ET}(X::DiscrGraph{ET}, β::Real, iters::Integer; seed = 167432777111, step::Integer = 1, hook = (x...)->true, C0::Union{Config,Void} = nothing,
                   staged_thr::Real = 0.5, staged_thr_fact::Real = 5.0)
    isfinite(β) || throw(ArgumentError("β must be finite, given: $β"))
    seed > 0 && srand(seed)
    Es = empty!(Array(ET, min(10^8, iters ÷ step)))

    N = getN(X)
    C::Config = C0 ≡ nothing ? Config(N) : C0
    C.N == N || throw(ArgumentError("Invalid C0, wrong N, expected $N, given: $(C.N)"))
    E = energy(X, C)
    ΔEcache = DeltaECache(X, C, β)
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
    println("samples = ", length(Es))
    println("iters = ", it)
    println("accept rate = ", accepted / it)
    println("frac. staged iters = ", staged_its / it)
    return Es, C
end

function rrrMC{ET}(X::DoubleGraph{ET}, β::Real, iters::Integer; seed = 167432777111, step::Integer = 1, hook = (x...)->true, C0::Union{Config,Void} = nothing,
                   staged_thr::Real = 0.5, staged_thr_fact::Real = 5.0)
    isfinite(β) || throw(ArgumentError("β must be finite, given: $β"))
    seed > 0 && srand(seed)
    Es = empty!(Array(ET, min(10^8, iters ÷ step)))

    N = getN(X)
    C::Config = C0 ≡ nothing ? Config(N) : C0
    C.N == N || throw(ArgumentError("Invalid C0, wrong N, expected $N, given: $(C.N)"))
    E = energy(X, C)
    X0 = discr_graph(X)
    ΔEcache = DeltaECache(X0, C, β)
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
    println("samples = ", length(Es))
    println("iters = ", it)
    println("accept rate = ", accepted / it)
    println("frac. staged iters = ", staged_its / it)
    return Es, C
end

### Begin BKL-related functions

function step_bkl(X::DiscrGraph, C::Config, ΔEcache::DeltaECache)
    skip = rand_skip(ΔEcache)
    move, ΔE = rand_move(ΔEcache)
    apply_move!(X, C, move, ΔEcache)
    #check_consistency(ΔEcache)
    return ΔE, skip
end

"""
    bklMC(X::DiscrGraph, β::Real, iters::Integer; keywords...)

Same as [`standardMC`](@ref), but uses the rejection-free method by Bortz, Kalos and Lebowitz. Each step takes more, but rejected moves
are essentially free, since they are skipped entirely.

The return values and the keyword arguments are the same as [`standardMC`](@ref), see the usage examples for that function. Note however
that this function can only be used with [`DiscrGraph`](@ref) models.

Note that the number of iterations includes the rejected moves. This makes the results directly comparable with those of `standardMC`. It also
means that increasing `β` at fixed `iters` will result in fewer steps being actually computed.
"""
function bklMC{ET}(X::DiscrGraph{ET}, β::Real, iters::Integer; seed = 167432777111, step::Integer = 1, hook = (x...)->true, C0::Union{Config,Void} = nothing)
    seed > 0 && srand(seed)
    Es = empty!(Array(ET, min(10^8, iters ÷ step)))

    N = getN(X)
    C::Config = C0 ≡ nothing ? Config(N) : C0
    C.N == N || throw(ArgumentError("Invalid C0, wrong N, expected $N, given: $(C.N)"))
    E = energy(X, C)
    ΔEcache = DeltaECache(X, C, β, false)
    #check_consistency(ΔEcache)

    it = 0
    accepted = 0
    nextstep = step
    stop = false
    while it < iters
        #@assert E == energy(X, C)
        #check_consistency(ΔEcache)

        ΔE, skip = step_bkl(X, C, ΔEcache)

        while it + skip + 1 ≥ nextstep
            push!(Es, E)
            hook(nextstep, X, C, accepted, E) || @goto out
            nextstep += step
            nextstep > iters && @goto out
        end
        it += skip + 1
        E += ΔE
        accepted += 1
    end
    @label out
    println("samples = ", length(Es))
    println("iters = ", it)
    println("accept rate = ", accepted / iters)
    println("true it = ", accepted)
    return Es, C
end

# function step_wtm!(X::AbstractGraph, C::Config, ΔEcache::THeap, β::Real, t0::Float64)
#     t, move = pick_next(ΔEcache)
#     δt = t - t0
#     ΔE = update_heap!(ΔEcache, X, C, β, move, t)
#     return ΔE, t, δt
# end

## function wtmMC{ET}(X::AbstractGraph{ET}, β::Real, iters::Integer; seed = 167432777111, step::Integer = 1, hook = (x...)->true, C0::Union{Config,Void} = nothing)
##     srand(seed)
##     Es = empty!(Array(ET, min(10^8, iters ÷ step)))
##
##     N = getN(X)
##     C::Config = C0 ≡ nothing ? Config(N) : C0
##     C.N == N || throw(ArgumentError("Invalid C0, wrong N, expected $N, given: $(C.N)"))
##     E = energy(X, C)
##     theap = THeap(X, C, β)
##     #check_consistency(ΔEcache)
##
##     it = 0
##     t = 0.0
##     nextstep = step
##     while it < iters
##         #@assert E == energy(X, C)
##
##         it += 1
##         t′, move = pick_next(theap)
##         δt = t′ - t
##         if (it % step == 0)
##             push!(Es, E)
##             hook(it, X, C, t, δt, E) || break
##         end
##         t = t′
##         ΔE = update_heap!(theap, X, C, β, move, t)
##         E += ΔE
##     end
##     println("global time = ", t)
##     println("ratio = ", t / it)
##     return Es, C
## end

# TODO: document!
function wtmMC{ET}(X::AbstractGraph{ET}, β::Real, samples::Integer; seed = 167432777111, step::Float64 = 1.0, hook = (x...)->true, C0::Union{Config,Void} = nothing)
    seed > 0 && srand(seed)
    Es = empty!(Array(ET, min(10^8, samples)))

    N = getN(X)
    C::Config = C0 ≡ nothing ? Config(N) : C0
    C.N == N || throw(ArgumentError("Invalid C0, wrong N, expected $N, given: $(C.N)"))
    E = energy(X, C)
    theap = THeap(X, C, β)
    #check_consistency(ΔEcache)

    num_moves = 0
    tmax = step * samples
    t = 0.0
    nextstep = step
    while t < tmax
        #@assert E == energy(X, C)

        t′, move = pick_next(theap)
        #δt = t′ - t
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
    println("samples = ", length(Es))
    println("num_moves = ", num_moves)
    println("global time = ", t)
    println("ratio = ", t / num_moves)
    return Es, C
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
        ΔEcache = DeltaECache(X, C, β)
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
        τ = Array(Float64, length(βr), 3)
        rr = Array(Float64, length(βr), 3)
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
