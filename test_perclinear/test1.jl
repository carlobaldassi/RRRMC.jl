module Test1

using RRRMC
include("Newton.jl")
using .Newton
#using Gaston
#set_terminal("qt")

using Interpolations

#include("aux_plot_stuff.jl")

function test1(;N = 1001, α = 0.6)
    P = round(Int, α * N)
    X = RRRMC.GraphPercLinear(N, P)

    optE = RRRMC.energy(X, optC) / N
    #optE = Float64(Float32(RRRMC.energy(X, optC)) / Float32(N))
    mag = 2 * (sum(optC.s) / N) - 1

    @show optE
    @show mag

    return X
end

function test_SA(;seed = 888, τ::Integer = 10^2, T0 = 3.0, T1 = 1e-15, C0 = nothing, N = 1001, α = 0.6)
    srand(seed + 96238575648268956374)
    P = round(Int, α * N)
    X = RRRMC.GraphPercLinear(N, P)

    tot_iters = τ * N

    iters = N
    num_steps = round(Int, tot_iters / iters)

    itst = round(Int, iters / 10)

    sseed = seed

    force = true
    function gen_hook(tag, β, fst, t0, it0)
        #isdir(dirname) || mkdir(dirname)
        #fn = gen_fname(tag, seed)
        fn = "test_SA_$(tag)_tau$(τ)_seed$(seed).txt"
        if fst
            !force && isfile(fn) && error("file $fn exists")
            f = open(fn, "w")
            println(f, "#it acc E time β")
            t0 = time()
        else
            isfile(fn) || error("file $fn not found")
            f = open(fn, "a")
        end
        hook = (it, X, C, acc, E) -> begin
            t = time() - t0
            #@assert abs(E - RRRMC.energy(X, C)) < 1e-10
            println(f, "$(it+it0) $acc $(E/N) $t $β")
            return true
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    fst = true; t0 = 0.0; C = C0; it0 = 0
    pp = randperm(N)
    #for β = linspace(β0, β1, num_steps)
    for T = linspace(T0, T1, num_steps)
        β = 1 / T
        info("β=$β")
        hook, cleanup, t0 = gen_hook("met", β, fst, t0, it0)
        try
            @time E, C = RRRMC.standardMC(X, β, iters, seed=sseed, step=itst, hook=hook, C0=C, pp=pp)
        finally
            cleanup()
        end
        fst = false
        it0 += iters
        #sseed += 624234
        #sseed += rand(1:10^9)
        sseed = 0
    end

    @show RRRMC.energy(X, C)

    return C, X
end

function read_equivβ_file(α::Float64, βQ::Float64, filename::AbstractString)
    βdict = Dict{Float64,Float64}()
    open(filename) do f
        for l in eachline(f)
            startswith(strip(l), "#") && continue
            ls = split(l)
            @assert length(ls) == 10
            parse(Float64, ls[1]) == α || continue
            parse(Float64, ls[2]) == βQ || continue
            Γ = parse(Float64, ls[3])
            β = parse(Float64, ls[5])
            haskey(βdict, Γ) && βdict[Γ] ≠ β && error("conflicting entry! α=$α βQ=$βQ Γ=$Γ: already read β=$(βdict[Γ]), now read β=$β")
            βdict[Γ] = β
        end
    end
    Γs = sort!(collect(keys(βdict)))
    βs = [βdict[Γ] for Γ in Γs]
    βitp = extrapolate(interpolate((Γs,), βs, Gridded(Linear())), Interpolations.Throw())
    return βitp
end

function test_SA_equivβ(;seed = 7001, τ::Integer = 10^2, C0 = nothing, N = 1001, α = 0.6, tag::AbstractString = "",
                         βQ::Float64 = 32.0, Γ0::Float64 = 2.5, Γ1::Float64 = 1e-2, num_steps = 30,
                         resfile = "equiv_beta_LinEnergy.txt")

    βitp = read_equivβ_file(α, βQ, resfile)

    srand(seed + 78123472837432554 + 96238575648268956374)
    P = round(Int, α * N)
    X = RRRMC.GraphPercLinear(N, P)

    tot_iters = τ * N
    iters = round(Int, tot_iters / num_steps)
    tot_iters % iters == 0 || warn("true tot_iters = $num_steps * $iters = $(num_steps * iters): effective τ=$(num_steps * iters / N) instead of $τ")

    itst = round(Int, iters / 100)

    sseed = seed

    force = true
    function gen_hook(alg, tag, β, Γ, fst, t0, it0)
        #isdir(dirname) || mkdir(dirname)
        #fn = gen_fname(alg, seed)
        fn = "test_SAequivbeta_$(alg)_tau$(τ)_seed$(seed)$(tag).txt"
        if fst
            !force && isfile(fn) && error("file $fn exists")
            f = open(fn, "w")
            println(f, "#it acc E time β (Γ)")
            t0 = time()
        else
            isfile(fn) || error("file $fn not found")
            f = open(fn, "a")
        end
        hook = (it, X, C, acc, E) -> begin
            t = time() - t0
            #@assert abs(E - RRRMC.energy(X, C)) < 1e-10
            println(f, "$(it+it0) $acc $(E/N) $t $β $Γ")
            return true
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    Γs = linspace(Γ0, Γ1, num_steps)
    βs = [βitp[Γ] for Γ in Γs]

    fst = true; t0 = 0.0; C = C0; it0 = 0
    # pp = randperm(N)
    for (β,Γ) in zip(βs, Γs)
        info("β=$β (Γ=$Γ)")
        islast = abs(Γ / Γ1 - 1) < 1e-10
        sc = islast ? 2 : 1
        hook, cleanup, t0 = gen_hook("met", tag, β, Γ, fst, t0, it0)
        try
            @time E, C = RRRMC.standardMC(X, β, iters ÷ sc, seed=sseed, step=itst, hook=hook, C0=C)
        finally
            cleanup()
        end
        if islast
            it0 += iters ÷ sc
            hook, cleanup, t0 = gen_hook("met", tag, βQ, 0.0, fst, t0, it0)
            try
                @time E, C = RRRMC.standardMC(X, β, iters - iters ÷ sc, seed=sseed, step=itst, hook=hook, C0=C)
            finally
                cleanup()
            end
            it0 += iters - iters ÷ sc
        else
            it0 += iters
        end
        fst = false
        #sseed += 624234
        #sseed += rand(1:10^9)
        sseed = 0
    end

    @show RRRMC.energy(X, C)

    return C, X
end

function C_from_C1(C1::RRRMC.Config, M::Integer)
    N = C1.N
    C = RRRMC.Config(N * M)
    for k = 1:M
        copy!(C.s, (k-1)*N + 1, C1.s, 1, N)
    end
    return C
end

function C1_from_C(C::RRRMC.Config, M::Integer)
    N = C.N
    @assert N % M == 0
    Nk = N ÷ M
    C1 = RRRMC.Config(Nk)
    avgm = 0.0
    unpol = 0
    for i = 1:Nk
        m = 0
        for k = 1:M
            m += 2 * C.s[(k-1)*Nk + i] - 1
        end
        avgm += abs(m / M)
        abs(m) ≠ M && (unpol += 1)
        C1.s[i] = m > 0 ? true : m < 0 ? false : rand(Bool)
    end
    avgm /= Nk
    frac_unpol = unpol / Nk
    @show avgm
    @show frac_unpol
    return C1
end

function test_QSA(;N::Integer = 1001, α::Real = 0.6, M::Integer = 32, alg = :rrr, τ::Integer = 10^2, seed = 78821000027346, tag::AbstractString = "", num_steps = 30,
                  β::Float64 = 32.0, β1::Float64 = NaN, precool::Bool = true, Γ0::Float64 = 2.5, Γ1::Float64 = 1e-2)
    isnan(β1) && (β1 = β)

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β, τ = 100, N = N, α = α, seed = seed + 78123472837432554)
    N = RRRMC.getN(Xs)
    P = Xs.P

    tot_iters = τ * N * M
    iters = round(Int, tot_iters / num_steps)
    tot_iters % iters == 0 || warn("true tot_iters = $num_steps * $iters = $(num_steps * iters): effective τ=$(num_steps * iters / (N * M)) instead of $τ")

    itst = round(Int, iters / 10)

    sseed = seed

    force = true
    function gen_hook(alg, tag, β, Γ, fst, t0, it0)
        #isdir(dirname) || mkdir(dirname)
        #fn = gen_fname(alg, seed)
        fn = "test_QSA_$(alg)_tau$(τ)_seed$(seed)$(tag).txt"
        if fst
            !force && isfile(fn) && error("file $fn exists")
            f = open(fn, "w")
            println(f, "#it acc E Ē time QE β Γ | q1s...")
            t0 = time()
        else
            isfile(fn) || error("file $fn not found")
            f = open(fn, "a")
        end
        hook = (it, X, C, acc, E) -> begin
            if isa(X, RRRMC.GraphQuant)
                Es = RRRMC.Renergies(X)
                mE = minimum(Es)
                aE = mean(Es)

                QE = RRRMC.Qenergy(X, C)
                q1 = RRRMC.overlaps(X)
            else
                mE = E
                aE = E
                QE = E
                q1 = ones(M÷2)
            end
            t = time() - t0
            println(f, "$(it+it0) $acc $(mE/N) $(aE/N) $t $QE $β $Γ | ", join(map(string, q1), " "))
            return true
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    fst = true; t0 = 0.0; C = nothing; it0 = 0
    precool && (C = C_from_C1(C1, M))
    local X
    for Γ = linspace(Γ0, Γ1, num_steps)
        info("Γ=$Γ β=$β")
        islast = abs(Γ / Γ1 - 1) < 1e-10
        sc = islast ? 2 : 1
        sMC = alg == :rrr ? RRRMC.rrrMC : RRRMC.standardMC
        hook, cleanup, t0 = gen_hook(string(alg), tag, β, Γ, fst, t0, it0)
        try
            X = RRRMC.GraphQPercT(Xs, M, Γ, β)
            @time E, C = sMC(X, β, iters ÷ sc, seed=sseed, step=itst, hook=hook, C0=C)
        finally
            cleanup()
        end
        if islast
            it0 += iters ÷ sc
            C1 = C1_from_C(C, M)
            hook, cleanup, t0 = gen_hook(string(alg), tag, β, 0.0, fst, t0, it0)
            try
                @time E, C1 = RRRMC.standardMC(Xs, β, iters - iters ÷ sc, seed=sseed, step=itst, hook=hook, C0=C1)
            finally
                cleanup()
            end
            C = C_from_C1(C1, M)
            it0 += iters - iters ÷ sc
        else
            it0 += iters
        end
        fst = false
        sseed = 0

        β += (β1 - β) * 0.1
    end

    Es = RRRMC.Renergies(X)
    @show Es
    @show mean(Es), std(Es)
    @show maximum(Es)
    @show minimum(Es)

    return C, Xs, X
end

function read_βΓ_file(α::Float64, filename::AbstractString)
    Γs = Float64[]
    βs = Float64[]
    open(filename) do f
        for l in eachline(f)
            startswith(strip(l), "#") && continue
            ls = split(l)
            @assert length(ls) == 12
            parse(Float64, ls[1]) == α || continue
            β = parse(Float64, ls[2])
            Γ = parse(Float64, ls[3])
            push!(βs, β)
            push!(Γs, Γ)
        end
    end
    @assert !isempty(length(Γs))
    p = sortperm(Γs)
    Γs = Γs[p]
    βs = βs[p]
    return Γs, βs
end

function getβfromΓ(Γ::Float64, Γs::Vector{Float64}, βs::Vector{Float64})
    L = length(Γs)
    @assert L > 0
    @assert length(βs) == L
    Γ0, β0 = Γs[1], βs[1]

    Γ ≤ Γ0 && return β0

    for i = 2:L
        Γ1, β1 = Γs[i], βs[i]
        Γ1 ≥ Γ && return β0 + ((Γ - Γ0) / (Γ1 - Γ0)) * (β1 - β0)
        Γ0, β0 = Γ1, β1
    end
    return β0
end

function test_QSA_fromfile(;N::Integer = 1001, α::Real = 0.6, M::Integer = 32, alg = :rrr, τ::Integer = 10^2, seed = 78821000027346, tag::AbstractString = "", num_steps = 30,
                            precool::Bool = false, Γ0::Float64 = 2.5, Γ1::Float64 = 1e-2, resfile::AbstractString = "results_QPerc_finiteT.fixed_sigma.txt")

    Γs, βs = read_βΓ_file(α, resfile)

    β = getβfromΓ(Γ0, Γs, βs)

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β, τ = 100, N = N, α = α, seed = seed + 78123472837432554)
    N = RRRMC.getN(Xs)
    P = Xs.P

    tot_iters = τ * N * M
    iters = round(Int, tot_iters / num_steps)
    tot_iters % iters == 0 || warn("true tot_iters = $num_steps * $iters = $(num_steps * iters): effective τ=$(num_steps * iters / (N * M)) instead of $τ")

    itst = round(Int, iters / 10)

    sseed = seed

    force = true
    function gen_hook(alg, tag, β, Γ, fst, t0, it0)
        fn = "test_QSA_fromfile_$(alg)_tau$(τ)_seed$(seed)$(tag).txt"
        if fst
            !force && isfile(fn) && error("file $fn exists")
            f = open(fn, "w")
            println(f, "#it acc E Ē time QE β Γ")
            t0 = time()
        else
            isfile(fn) || error("file $fn not found")
            f = open(fn, "a")
        end
        hook = (it, X, C, acc, E) -> begin
            if isa(X, RRRMC.GraphQuant)
                Es = RRRMC.Renergies(X)
                mE = minimum(Es)
                aE = mean(Es)

                QE = RRRMC.Qenergy(X, C)
                q1 = RRRMC.overlaps(X)
            else
                mE = E
                aE = E
                QE = E
                q1 = ones(M÷2)
            end
            t = time() - t0
            println(f, "$(it+it0) $acc $(mE/N) $(aE/N) $t $QE $β $Γ | ", join(map(string, q1), " "))
            return true
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    fst = true; t0 = 0.0; C = nothing; it0 = 0
    precool && (C = C_from_C1(C1, M))
    local X
    for Γ = linspace(Γ0, Γ1, num_steps)
        β = getβfromΓ(Γ, Γs, βs)
        info("Γ=$Γ β=$β")
        islast = abs(Γ / Γ1 - 1) < 1e-10
        sc = islast ? 2 : 1
        sMC = alg == :rrr ? RRRMC.rrrMC : RRRMC.standardMC
        hook, cleanup, t0 = gen_hook(string(alg), tag, β, Γ, fst, t0, it0)
        try
            X = RRRMC.GraphQPercT(Xs, M, Γ, β)
            @time E, C = sMC(X, β, iters ÷ sc, seed=sseed, step=itst, hook=hook, C0=C)
        finally
            cleanup()
        end
        if islast
            it0 += iters ÷ sc
            C1 = C1_from_C(C, M)
            hook, cleanup, t0 = gen_hook(string(alg), tag, β, 0.0, fst, t0, it0)
            try
                @time E, C1 = RRRMC.standardMC(Xs, β, iters - iters ÷ sc, seed=sseed, step=itst, hook=hook, C0=C1)
            finally
                cleanup()
            end
            C = C_from_C1(C1, M)
            it0 += iters - iters ÷ sc
        else
            it0 += iters
        end
        fst = false
        sseed = 0
    end

    Es = RRRMC.Renergies(X)
    @show Es
    @show mean(Es), std(Es)
    @show maximum(Es)
    @show minimum(Es)

    return C, Xs, X
end

# Robust Ensemble

function C_from_C1_tr(C1::RRRMC.Config, M::Integer)
    Nk = C1.N
    C = RRRMC.Config(Nk * M)
    s = C.s
    s1 = C1.s
    for k = 1:M
        for (i,j) = enumerate(k:M:(k + M * (Nk-1)))
            s[j] = s1[i]
        end
    end
    return C
end

function C1_from_C_tr(C::RRRMC.Config, M::Integer)
    N = C.N
    @assert N % M == 0
    Nk = N ÷ M
    C1 = RRRMC.Config(Nk)
    avgm = 0.0
    unpol = 0
    for i = 1:Nk
        m = 0
        for k = 1:M
            m += 2 * C.s[(i-1)*M + k] - 1
        end
        avgm += abs(m / M)
        abs(m) ≠ M && (unpol += 1)
        C1.s[i] = m > 0 ? true : m < 0 ? false : rand(Bool)
    end
    avgm /= Nk
    frac_unpol = unpol / Nk
    @show avgm
    @show frac_unpol
    return C1
end

function C1_from_C_tr(C::RRRMC.Config, M::Integer, X::RRRMC.GraphRobustEnsemble)
    N = C.N
    @assert N % M == 0
    Nk = N ÷ M
    C1 = RRRMC.Config(Nk)
    Es = RRRMC.REenergies(X)

    _, k = findmin(Es)
    s = C.s
    s1 = C1.s
    for (i,j) = enumerate(k:M:(k + M * (Nk-1)))
        s1[i] = s[j]
    end
    return C1
end

function test_RSA(;N::Integer = 1001, α::Real = 0.6, M::Integer = 32, alg = :rrr, τ::Integer = 10^2, seed = 78821000027346, tag::AbstractString = "", num_steps::Integer = 30,
                  β::Float64 = 1.0, β1::Float64 = 0.1, βsched = :additive, precool::Bool = true,
                  γ0::Float64 = 0.1, γ1::Float64 = 1e5)

    @assert βsched ∈ [:additive, :multiplicative]

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β, τ = 100, N = N, α = α, seed = seed + 78123472837432554)
    N = RRRMC.getN(Xs)
    P = Xs.P

    tot_iters = τ * N
    iters = round(Int, tot_iters / num_steps)
    tot_iters % iters == 0 || warn("true tot_iters = $num_steps * $iters = $(num_steps * iters): effective τ=$(num_steps * iters / N)")

    itst = round(Int, iters / 10)

    if βsched == :additive
        βs = (β1 - β) / (num_steps - 1)
    else # βsched == :multiplicative
        βf = (β1/β)^(1/(num_steps-1))
    end

    g0 = -Float64(log(tanh(M/β1 * big(γ0))))
    g1 = -Float64(log(tanh(M/β1 * big(γ1))))

    sseed = seed

    force = true
    function gen_hook(alg, tag, β, γ, fst, t0, it0)
        #isdir(dirname) || mkdir(dirname)
        fn = "test_RSA_$(alg)_tau$(τ)_seed$(seed)$(tag).txt"
        if fst
            !force && isfile(fn) && error("file $fn exists")
            f = open(fn, "w")
            println(f, "#it acc Emin Emax <E> σE time β γ")
            t0 = time()
        else
            isfile(fn) || error("file $fn not found")
            f = open(fn, "a")
        end
        hook = (it, X, C, acc, E) -> begin
            if isa(X, RRRMC.GraphRobustEnsemble)
                Es = RRRMC.REenergies(X)
                mE, ME = extrema(Es)
                aE, sE = mean(Es), std(Es)
            else
                mE, ME, aE, sE = E, E, E, 0.0
            end
            t = time() - t0
            println(f, "$(it+it0) $acc $(mE/N) $(ME/N) $(aE/N) $(sE/N) $t $β $γ")
            return true
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    userrr = true

    fst = true; t0 = 0.0; C = nothing; it0 = 0
    precool && (C = C_from_C1_tr(C1, M))
    local X
    for g = linspace(g0, g1, num_steps)
        #γ = exp(g)
        γ = min(γ1, Float64(atanh(exp(-big(g))) * β1/M))
        info("g=$(g) γ=$(γ) β=$β useRRR=$userrr")
        islast = abs(g / g1 - 1) < 1e-10
        sc = islast ? 2 : 1
        sMC = alg == :rrr ? RRRMC.rrrMC : RRRMC.standardMC
        hook, cleanup, t0 = gen_hook(string(alg), tag, β, γ, fst, t0, it0)
        try
            X = RRRMC.GraphPercRE(Xs, M, γ, β)
            @time E, C = sMC(X, β, iters ÷ sc, seed=sseed, step=itst, hook=hook, C0=C)
        finally
            cleanup()
        end
        if islast
            it0 += iters ÷ sc
            C1 = C1_from_C_tr(C, M, X)
            hook, cleanup, t0 = gen_hook(string(alg), tag, 20.0, 0.0, fst, t0, it0)
            try
                @time E, C1 = RRRMC.standardMC(Xs, 20.0, iters - iters ÷ sc, seed=sseed, step=itst, hook=hook, C0=C1)
            finally
                cleanup()
            end
            C = C_from_C1_tr(C1, M)
            it0 += iters - iters ÷ sc
        else
            it0 += iters
        end
        fst = false
        sseed = 0

        if βsched == :additive
            β += βs
        else # βsched == :multiplicative
            β *= βf
        end
    end

    Es = RRRMC.REenergies(X)
    @show Es
    @show mean(Es), std(Es)
    @show maximum(Es)
    @show minimum(Es)
    @show RRRMC.energy(Xs, C1)

    return C, Xs, X
end

# Local Entropy

function C_from_C1_tr1(C1::RRRMC.Config, M::Integer)
    Nk = C1.N
    C = RRRMC.Config(Nk * (M+1))
    s = C.s
    s1 = C1.s
    for k = 1:(M+1)
        for (i,j) = enumerate(k:(M+1):(k + (M+1) * (Nk-1)))
            s[j] = s1[i]
        end
    end
    return C
end

function C1_from_C_tr1(C::RRRMC.Config, M::Integer)
    N = C.N
    @assert N % (M+1) == 0
    Nk = N ÷ (M+1)
    C1 = RRRMC.Config(Nk)
    s = C.s
    s1 = C1.s
    for (i,j) = enumerate(1:(M+1):(1 + (M+1) * (Nk-1)))
        s1[i] = s[j]
    end
    return C1
end


function test_LSA(;N::Integer = 1001, α::Real = 0.6, M::Integer = 32, alg = :rrr, τ::Integer = 10^2, seed = 78821000027346, tag::AbstractString = "", num_steps::Integer = 30,
                  β::Float64 = 0.1, β1::Float64 = 1.5, precool::Bool = false, init_equal = false,
                  γ0::Float64 = 0.1, γ1::Float64 = 1e5)

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β1, τ = precool ? 100 : 2, N = N, α = α, seed = seed + 78123472837432554)
    N = RRRMC.getN(Xs)
    P = Xs.P

    tot_iters = τ * N
    iters = round(Int, tot_iters / num_steps)
    tot_iters % iters == 0 || warn("true tot_iters = $num_steps * $iters = $(num_steps * iters): effective τ=$(num_steps * iters / N)")

    itst = round(Int, iters / 10)

    βf = (β1/β)^(1/(num_steps-1))

    #β1 = β + num_steps * β

    #g0 = -Float64(log(tanh(M/β1 * big(γ0))))
    #g1 = -Float64(log(tanh(M/β1 * big(γ1))))
    g0 = -Float64(log(tanh(M * big(γ0))))
    g1 = -Float64(log(tanh(M * big(γ1))))
    #g0 = log(γ0)
    #g1 = log(γ1)

    sseed = seed

    force = true
    function gen_hook(alg, tag, β, γ, fst, t0, it0)
        #isdir(dirname) || mkdir(dirname)
        fn = "test_LSA_$(alg)_tau$(τ)_seed$(seed)$(tag).txt"
        if fst
            !force && isfile(fn) && error("file $fn exists")
            f = open(fn, "w")
            println(f, "#it acc Emin Emax <E> σE time β γ")
            t0 = time()
        else
            isfile(fn) || error("file $fn not found")
            f = open(fn, "a")
        end
        hook = (it, X, C, acc, E) -> begin
            if isa(X, RRRMC.GraphLocEntr)
                Es = RRRMC.LEenergies(X)
                cE = RRRMC.cenergy(X)
                mE, ME = extrema(Es)
                aE, sE = mean(Es), std(Es)
            else
                mE, ME, aE, sE, cE = E, E, E, 0.0, E
            end
            t = time() - t0
            println(f, "$(it+it0) $acc $(mE/N) $(ME/N) $(aE/N) $(sE/N) $(cE/N) $t $β $γ")
            return true
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    userrr = true

    fst = true; t0 = 0.0; C = nothing; it0 = 0
    !precool && init_equal && (C1 = RRRMC.Config(N))
    (precool || init_equal) && (C = C_from_C1_tr1(C1, M))
    local X
    for g = linspace(g0, g1, num_steps)
        #γ = exp(g)
        #γ = min(γ1, Float64(atanh(exp(-big(g))) * β1/M))
        γ = min(γ1, Float64(atanh(exp(-big(g))) / M))
        info("g=$(g) γ=$(γ) β=$β useRRR=$userrr")
        islast = abs(g / g1 - 1) < 1e-10
        sc = islast ? 2 : 1
        sMC = alg == :rrr ? RRRMC.rrrMC : RRRMC.standardMC
        hook, cleanup, t0 = gen_hook(string(alg), tag, β, γ, fst, t0, it0)
        try
            X = RRRMC.GraphPercLE(Xs, M, γ, β)
            @time E, C = sMC(X, β, iters ÷ sc, seed=sseed, step=itst, hook=hook, C0=C)
        finally
            cleanup()
        end
        if islast
            it0 += iters ÷ sc
            C1 = C1_from_C_tr1(C, M)
            hook, cleanup, t0 = gen_hook(string(alg), tag, 20.0, 0.0, fst, t0, it0)
            try
                @time E, C1 = RRRMC.standardMC(Xs, 20.0, iters - iters ÷ sc, seed=sseed, step=itst, hook=hook, C0=C1)
            finally
                cleanup()
            end
            C = C_from_C1_tr1(C1, M)
            it0 += iters - iters ÷ sc
        else
            it0 += iters
        end
        fst = false
        sseed = 0

        #β += (β1 - β) * βf
        #β += βs
        β *= βf
    end

    Es = RRRMC.LEenergies(X)
    @show Es
    @show mean(Es), std(Es)
    @show maximum(Es)
    @show minimum(Es)
    @show RRRMC.energy(Xs, C1)

    return C, Xs, X
end

end # module
