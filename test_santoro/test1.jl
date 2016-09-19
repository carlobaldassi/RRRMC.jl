module Test1

using RRRMC
include("Newton.jl")
using .Newton
#using Gaston
#set_terminal("qt")

#include("aux_plot_stuff.jl")

function test1()
    #fname = "80x80_couplings.dat"
    #fname = "80x80_uniform1.txt"
    fname = "1000_uniform_systems/text_files_and_gs/80x80_uniform1.txt"
    X = RRRMC.GraphEAContSimple(fname)
    N = RRRMC.getN(X)
    optC = RRRMC.Config(N)
    fill!(optC.s, 0)

    open("conf_uniform1.txt") do f
        for l in eachline(f)
            i = parse(Int, l)
            optC.s[i] = 1
        end
    end

    #flipbits!(optC.s)

    optE = RRRMC.energy(X, optC) / N
    #optE = Float64(Float32(RRRMC.energy(X, optC)) / Float32(N))
    mag = 2 * (sum(optC.s) / N) - 1

    @show optE
    @show mag

    return X
end

function test_SA(;seed = 888, τ::Integer = 10^2, T0 = 3.0, T1 = 1e-15, C0 = nothing)
    #fname = "80x80_couplings.dat"
    fname = "1000_uniform_systems/text_files_and_gs/80x80_uniform1.txt"
    X = RRRMC.GraphEAContSimple(fname)
    N = RRRMC.getN(X)

    #optE = -1.5805166767932775
    #optE = -1.59251261170312 #???
    optE = -1.5925126871330073
    #optE = -1.593443042567935

    #T0 = 3.0
    #T1 = 1e-15

    #β0 = 1/3
    #β1 = 10.0

    tot_iters = τ * N

    iters = N
    num_steps = round(Int, tot_iters / iters)

    itst = round(Int, iters / 10)

    #srand(sseed)
    #C = RRRMC.Config(RRRMC.getN(X))
    #allE = Float64[]

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

function test_QSA(;M::Integer = 32, alg = :rrr, τ::Integer = 10^2, seed = 78821000027346, tag::AbstractString = "", num_steps = 30,
                  β::Float64 = 32.0, β1::Float64 = NaN, precool::Bool = true, Γ0::Float64 = 2.5, Γ1::Float64 = 1e-2)
    #fname = "80x80_couplings.dat"
    #fname = "1000_uniform_systems/text_files_and_gs/80x80_uniform1.txt"
    #Xs = RRRMC.GraphEAContSimple(fname)
    #N = RRRMC.getN(Xs)

    #optE = -1.5805166767932775
    #optE = -1.5925126871330073

    isnan(β1) && (β1 = β)

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β, τ = 100)
    N = RRRMC.getN(Xs)

    tot_iters = τ * N
    iters = round(Int, tot_iters / num_steps)
    tot_iters % iters == 0 || warn("true tot_iters = $num_steps * $iters = $(num_steps * iters): effective τ=$(num_steps * iters / N)")

    itst = round(Int, iters / 10)

    #@assert τ ≥ 100
    #λ = (1 - 100/τ)^10 # very empirical formula...
    #accrthresh = 0.07

    sseed = seed

    #accrate = Float64[0.1]
    force = true
    function gen_hook(alg, tag, β, Γ, fst, t0, it0)
        #isdir(dirname) || mkdir(dirname)
        #fn = gen_fname(alg, seed)
        fn = "test_QSA_$(alg)_tau$(τ)_seed$(seed)$(tag).txt"
        if fst
            !force && isfile(fn) && error("file $fn exists")
            f = open(fn, "w")
            println(f, "#it acc E time QE β Γ")
            t0 = time()
        else
            isfile(fn) || error("file $fn not found")
            f = open(fn, "a")
        end
        hook = (it, X, C, acc, E) -> begin
            if isa(X, RRRMC.GraphQuant)
                mE = minimum(RRRMC.Renergies(X))
                QE = RRRMC.Qenergy(X, C)
            else
                mE = E
                QE = E
            end
            t = time() - t0
            println(f, "$(it+it0) $acc $(mE/N) $t $QE $β $Γ")
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
            X = RRRMC.GraphQEAT(Xs, M, Γ, β)
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
    @show Es/N - optE
    @show mean(Es/N - optE), std(Es/N)
    @show maximum(Es)/N - optE
    @show minimum(Es)/N - optE
    return C, Xs, X
end

findb(b, it, ns, x) = abs(b * (-1 + (b / (b + it))^(-1/ns)) - x)

function test_QSA_alt(;M::Integer = 32, alg = :rrr, τ::Integer = 10^2, seed = 78821000027346, tag::AbstractString = "", num_steps::Integer = 30,
                      β::Float64 = 32.0, β1::Float64 = NaN, precool::Bool = true,
                      Γ0::Float64 = 10.0, Γ1::Float64 = 1e-2)
    #fname = "80x80_couplings.dat"
    #fname = "1000_uniform_systems/text_files_and_gs/80x80_uniform1.txt"
    #Xs = RRRMC.GraphEAContSimple(fname)
    #N = RRRMC.getN(Xs)

    #optE = -1.5805166767932775
    #optE = -1.5925126871330073

    isnan(β1) && (β1 = β)

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β, τ = 100)
    N = RRRMC.getN(Xs)

    tot_iters = τ * N

    nres = newton(b->findb(b, tot_iters, num_steps, N*M), Float64(N*M), NewtonParameters(1e-3, 1e-3, 0, 10^8))
    @assert nres[1] == true
    b = nres[2]
    #error()
    @show b
    iters = diff(map(x->round(Int,exp(x)), linspace(log(b), log(tot_iters+b), num_steps+1)))
    @show iters
    @show sum(iters)

    # x = N * M
    # nres = newton(ρ->abs((ρ^num_steps - 1)/(ρ - 1) - tot_iters/x), 1.1, NewtonParameters(1e-3, 1e-3, 0, 10^8))
    # @assert nres[1] == true
    # ρ = nres[2]
    # @show ρ
    # iters = [round(Int, x * ρ^t) for t = 0:(num_steps-1)]
    # @show iters
    # @show sum(iters)

    tot_iters ≠ sum(iters) && warn("true tot_iters = $(sum(iters)): effective τ=$(sum(iters) / N)")

    itst = map(it->max(1,round(Int, it / 10)), iters)
    @show itst
    Γ = linspace(Γ0, Γ1, num_steps)
    #Γ = collect(linspace(Γ0, Γ1, num_steps+1))[1:end-1]

    #@assert τ ≥ 100
    #λ = (1 - 100/τ)^10 # very empirical formula...
    #accrthresh = 1.0 # 0.07

    sseed = seed

    force = true
    function gen_hook(alg, tag, β, Γ, fst, t0, it0)
        #isdir(dirname) || mkdir(dirname)
        fn = "test_QSA_$(alg)_tau$(τ)_seed$(seed).ALT$(tag).txt"
        if fst
            !force && isfile(fn) && error("file $fn exists")
            f = open(fn, "w")
            println(f, "#it acc E time QE β Γ")
            t0 = time()
        else
            isfile(fn) || error("file $fn not found")
            f = open(fn, "a")
        end
        hook = (it, X, C, acc, E) -> begin
            mE = minimum(RRRMC.Renergies(X))
            QE = RRRMC.Qenergy(X, C)
            t = time() - t0
            println(f, "$(it+it0) $acc $(mE/N) $t $QE $β $Γ")
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    fst = true; t0 = 0.0; C = nothing; it0 = 0
    precool && (C = C_from_C1(C1, M))
    local X
    for s = 1:num_steps
        info("Γ=$(Γ[s]) β=$β")
        sMC = (alg == :rrr && userrr) ? RRRMC.rrrMC : RRRMC.standardMC
        hook, cleanup, t0 = gen_hook(string(alg), tag, β, Γ[s], fst, t0, it0)
        try
            X = RRRMC.GraphQEAT(Xs, M, Γ[s], β)
            #@show X, β, iters[s], sseed, itst[s]
            @time E, C = sMC(X, β, iters[s], seed=sseed, step=itst[s], hook=hook, C0=C)
        finally
            cleanup()
        end
        fst = false
        it0 += iters[s]
        #sseed += 624234
        sseed = 0

        β += (β1 - β) * 0.1
    end

    Es = RRRMC.Renergies(X)
    @show Es/N - optE
    @show mean(Es/N - optE), std(Es/N)
    @show maximum(Es)/N - optE
    @show minimum(Es)/N - optE

    return C, Xs, X
end

function test_QSA_alt2(;M::Integer = 32, alg = :rrr, τ::Integer = 10^2, seed = 78821000027346, tag::AbstractString = "", num_steps::Integer = 30,
                       β::Float64 = 32.0, β1::Float64 = NaN, precool::Bool = true, Γ0::Float64 = 2.5, Γ1::Float64 = 1e-2)
    isnan(β1) && (β1 = β)

    #fname = "1000_uniform_systems/text_files_and_gs/80x80_uniform1.txt"
    optE = -1.5925126871330073

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β, τ = 100)
    N = RRRMC.getN(Xs)

    tot_iters = τ * N
    iters = round(Int, tot_iters / num_steps)
    tot_iters % iters == 0 || warn("true tot_iters = $num_steps * $iters = $(num_steps * iters): effective τ=$(num_steps * iters / N)")

    itst = round(Int, iters / 10)

    g0 = log(1/2β * log(coth(β * Γ0 / M)))
    g1 = log(1/2β * log(coth(β * Γ1 / M)))

    sseed = seed

    force = true
    function gen_hook(alg, tag, β, Γ, fst, t0, it0)
        #isdir(dirname) || mkdir(dirname)
        fn = "test_QSA_$(alg)_tau$(τ)_seed$(seed).ALT2$(tag).txt"
        if fst
            !force && isfile(fn) && error("file $fn exists")
            f = open(fn, "w")
            println(f, "#it acc E time QE β Γ")
            t0 = time()
        else
            isfile(fn) || error("file $fn not found")
            f = open(fn, "a")
        end
        hook = (it, X, C, acc, E) -> begin
            if isa(X, RRRMC.GraphQuant)
                mE = minimum(RRRMC.Renergies(X))
                QE = RRRMC.Qenergy(X, C)
            else
                mE = E
                QE = E
            end
            t = time() - t0
            println(f, "$(it+it0) $acc $(mE/N) $t $QE $β $Γ")
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    userrr = true

    fst = true; t0 = 0.0; C = nothing; it0 = 0
    precool && (C = C_from_C1(C1, M))
    local X
    for g = linspace(g0, g1, num_steps)
        γ = exp(g)
        Γ = Float64(M / β * log(coth(0.5*acosh(exp(2β * big(γ))))))
        info("Γ=$(Γ) g=$(g) 4γ=$(4γ) β=$β useRRR=$userrr")
        islast = abs(g / g1 - 1) < 1e-10
        sc = islast ? 2 : 1
        sMC = alg == :rrr ? RRRMC.rrrMC : RRRMC.standardMC
        hook, cleanup, t0 = gen_hook(string(alg), tag, β, Γ, fst, t0, it0)
        try
            X = RRRMC.GraphQEAT(Xs, M, Γ, β)
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
    @show Es/N - optE
    @show mean(Es/N - optE), std(Es/N)
    @show maximum(Es)/N - optE
    @show minimum(Es)/N - optE

    return C, Xs, X
end

function test_QSA_alt3(;M::Integer = 32, alg = :rrr, τ::Integer = 10^2, seed = 78821000027346, tag::AbstractString = "", num_steps::Integer = 30,
                       β::Float64 = 32.0, β1::Float64 = NaN, precool::Bool = true, Γ0::Float64 = 2.5, Γ1::Float64 = 1e-2)
    isnan(β1) && (β1 = β)

    #fname = "1000_uniform_systems/text_files_and_gs/80x80_uniform1.txt"
    optE = -1.5925126871330073

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β, τ = 100)
    N = RRRMC.getN(Xs)

    tot_iters = τ * N
    iters = round(Int, tot_iters / num_steps)
    tot_iters % iters == 0 || warn("true tot_iters = $num_steps * $iters = $(num_steps * iters): effective τ=$(num_steps * iters / N)")

    itst = round(Int, iters / 10)

    γ0 = 1/2β * log(coth(β * Γ0 / M))
    γ1 = 1/2β * log(coth(β * Γ1 / M))

    sseed = seed

    force = true
    function gen_hook(alg, tag, β, Γ, fst, t0, it0)
        #isdir(dirname) || mkdir(dirname)
        fn = "test_QSA_$(alg)_tau$(τ)_seed$(seed).ALT3$(tag).txt"
        if fst
            !force && isfile(fn) && error("file $fn exists")
            f = open(fn, "w")
            println(f, "#it acc E time QE β Γ")
            t0 = time()
        else
            isfile(fn) || error("file $fn not found")
            f = open(fn, "a")
        end
        hook = (it, X, C, acc, E) -> begin
            if isa(X, RRRMC.GraphQuant)
                mE = minimum(RRRMC.Renergies(X))
                QE = RRRMC.Qenergy(X, C)
            else
                mE = E
                QE = E
            end
            t = time() - t0
            println(f, "$(it+it0) $acc $(mE/N) $t $QE $β $Γ")
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    userrr = true

    fst = true; t0 = 0.0; C = nothing; it0 = 0
    precool && (C = C_from_C1(C1, M))
    local X
    for γ = linspace(γ0, γ1, num_steps)
        Γ = Float64(M / β * log(coth(0.5*acosh(exp(2β * big(γ))))))
        info("Γ=$(Γ) 4γ=$(4γ) β=$β useRRR=$userrr")
        islast = abs(γ / γ1 - 1) < 1e-10
        sc = islast ? 2 : 1
        sMC = alg == :rrr ? RRRMC.rrrMC : RRRMC.standardMC
        hook, cleanup, t0 = gen_hook(string(alg), tag, β, Γ, fst, t0, it0)
        try
            X = RRRMC.GraphQEAT(Xs, M, Γ, β)
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
    @show Es/N - optE
    @show mean(Es/N - optE), std(Es/N)
    @show maximum(Es)/N - optE
    @show minimum(Es)/N - optE

    return C, Xs, X
end

function test_QSA_alt4(;M::Integer = 32, alg = :rrr, τ::Integer = 10^2, seed = 78821000027346, tag::AbstractString = "", num_steps::Integer = 30,
                       β::Float64 = 32.0, β1::Float64 = NaN, precool::Bool = true, Γ0::Float64 = 2.5, Γ1::Float64 = 1e-2)
    isnan(β1) && (β1 = β)

    #fname = "1000_uniform_systems/text_files_and_gs/80x80_uniform1.txt"
    optE = -1.5925126871330073

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β, τ = 100)
    N = RRRMC.getN(Xs)

    tot_iters = τ * N
    iters = round(Int, tot_iters / num_steps)
    tot_iters % iters == 0 || warn("true tot_iters = $num_steps * $iters = $(num_steps * iters): effective τ=$(num_steps * iters / N)")

    itst = round(Int, iters / 10)

    g0 = Γ0^2
    g1 = Γ1^2

    sseed = seed

    force = true
    function gen_hook(alg, tag, β, Γ, fst, t0, it0)
        #isdir(dirname) || mkdir(dirname)
        fn = "test_QSA_$(alg)_tau$(τ)_seed$(seed).ALT4$(tag).txt"
        if fst
            !force && isfile(fn) && error("file $fn exists")
            f = open(fn, "w")
            println(f, "#it acc E time QE β Γ")
            t0 = time()
        else
            isfile(fn) || error("file $fn not found")
            f = open(fn, "a")
        end
        hook = (it, X, C, acc, E) -> begin
            if isa(X, RRRMC.GraphQuant)
                mE = minimum(RRRMC.Renergies(X))
                QE = RRRMC.Qenergy(X, C)
            else
                mE = E
                QE = E
            end
            t = time() - t0
            println(f, "$(it+it0) $acc $(mE/N) $t $QE $β $Γ")
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    userrr = true

    fst = true; t0 = 0.0; C = nothing; it0 = 0
    precool && (C = C_from_C1(C1, M))
    local X
    for g = linspace(g0, g1, num_steps)
        #γ = (g)^10
        #Γ = Float64(M / β * log(coth(0.5*acosh(exp(2β * big(γ))))))
        Γ = g^0.5
        info("Γ=$(Γ) g=$(g) 4γ=$(4γ) β=$β useRRR=$userrr")
        islast = abs(g / g1 - 1) < 1e-10
        sc = islast ? 2 : 1
        sMC = alg == :rrr ? RRRMC.rrrMC : RRRMC.standardMC
        hook, cleanup, t0 = gen_hook(string(alg), tag, β, Γ, fst, t0, it0)
        try
            X = RRRMC.GraphQEAT(Xs, M, Γ, β)
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
    @show Es/N - optE
    @show mean(Es/N - optE), std(Es/N)
    @show maximum(Es)/N - optE
    @show minimum(Es)/N - optE

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

function test_RSA(;M::Integer = 32, alg = :rrr, τ::Integer = 10^2, seed = 78821000027346, tag::AbstractString = "", num_steps::Integer = 30,
                  β::Float64 = 1.0, βs::Float64 = 0.1, precool::Bool = true,
                  γ0::Float64 = 0.1, γ1::Float64 = 1e5)

    #fname = "1000_uniform_systems/text_files_and_gs/80x80_uniform1.txt"
    optE = -1.5925126871330073

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β, τ = precool ? 100 : 2)
    N = RRRMC.getN(Xs)

    tot_iters = τ * N
    iters = round(Int, tot_iters / num_steps)
    tot_iters % iters == 0 || warn("true tot_iters = $num_steps * $iters = $(num_steps * iters): effective τ=$(num_steps * iters / N)")

    itst = round(Int, iters / 10)

    β1 = β + num_steps * β

    #g0 = log(γ0)
    #g1 = log(γ1)
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
            if isa(X, RRRMC.GraphRepl)
                Es = RRRMC.REenergies(X)
                mE, ME = extrema(Es)
                aE, sE = mean(Es), std(Es)
            else
                mE, ME, aE, sE = E, E, E, Inf
            end
            t = time() - t0
            println(f, "$(it+it0) $acc $(mE/N) $(ME/N) $(aE/N) $(sE/N) $t $β $γ")
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
            X = RRRMC.GraphEARE(Xs, M, γ, β)
            @time E, C = sMC(X, β, iters ÷ sc, seed=sseed, step=itst, hook=hook, C0=C)
        finally
            cleanup()
        end
        if islast
            it0 += iters ÷ sc
            C1 = C1_from_C_tr(C, M)
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

        #β += (β1 - β) * βf
        β += βs
    end

    Es = RRRMC.REenergies(X)
    @show Es/N - optE
    @show mean(Es/N - optE), std(Es/N)
    @show maximum(Es)/N - optE
    @show minimum(Es)/N - optE

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


function test_LSA(;M::Integer = 32, alg = :rrr, τ::Integer = 10^2, seed = 78821000027346, tag::AbstractString = "", num_steps::Integer = 30,
                  β::Float64 = 1.0, βs::Float64 = 0.1, precool::Bool = true,
                  γ0::Float64 = 0.1, γ1::Float64 = 1e5)

    #fname = "1000_uniform_systems/text_files_and_gs/80x80_uniform1.txt"
    optE = -1.5925126871330073

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β, τ = precool ? 100 : 2)
    N = RRRMC.getN(Xs)

    tot_iters = τ * N
    iters = round(Int, tot_iters / num_steps)
    tot_iters % iters == 0 || warn("true tot_iters = $num_steps * $iters = $(num_steps * iters): effective τ=$(num_steps * iters / N)")

    itst = round(Int, iters / 10)

    β1 = β + num_steps * β

    #g0 = log(γ0)
    #g1 = log(γ1)
    g0 = -Float64(log(tanh(M/β1 * big(γ0))))
    g1 = -Float64(log(tanh(M/β1 * big(γ1))))

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
                mE, ME = extrema(Es)
                aE, sE = mean(Es), std(Es)
            else
                mE, ME, aE, sE = E, E, E, Inf
            end
            t = time() - t0
            println(f, "$(it+it0) $acc $(mE/N) $(ME/N) $(aE/N) $(sE/N) $t $β $γ")
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    userrr = true

    fst = true; t0 = 0.0; C = nothing; it0 = 0
    precool && (C = C_from_C1_tr1(C1, M))
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
            X = RRRMC.GraphEALE(Xs, M, γ, β)
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
        β += βs
    end

    Es = RRRMC.LEenergies(X)
    @show Es/N - optE
    @show mean(Es/N - optE), std(Es/N)
    @show maximum(Es)/N - optE
    @show minimum(Es)/N - optE

    return C, Xs, X
end

end # module
