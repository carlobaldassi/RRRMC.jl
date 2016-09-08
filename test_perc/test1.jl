module Test1

using RRRMC
include("Newton.jl")
using .Newton
#using Gaston
#set_terminal("qt")

#include("aux_plot_stuff.jl")

function test1(;N = 1001, α = 0.6)
    P = round(Int, α * N)
    X = RRRMC.GraphPerc(N, P)

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
    X = RRRMC.GraphPerc(N, P)

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

function test_QSA(;N::Integer = 1001, α::Real = 0.6, M::Integer = 32, tag = :rrr, τ::Integer = 10^2, seed = 78821000027346)
    β = 32.0
    β1 = 32.0

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β, τ = 100, N = N, α = α, seed = seed + 78123472837432554)
    N = RRRMC.getN(Xs)
    P = Xs.P

    Γ0 = 2.5
    Γ1 = 1e-2

    tot_iters = τ * N
    iters = N * M * max(1, τ ÷ 1000)


    num_steps = round(Int, tot_iters / iters)
    tot_iters % iters == 0 || warn("true tot_iters = $num_steps * $iters = $(num_steps * iters): effective τ=$(num_steps * iters / N)")

    itst = round(Int, iters / 10)

    @assert τ ≥ 100
    λ = (1 - 100/τ)^10 # very empirical formula...
    accrthresh = 0.07

    sseed = seed

    accrate = Float64[0.1]
    force = true
    function gen_hook(tag, β, Γ, fst, t0, it0)
        #isdir(dirname) || mkdir(dirname)
        #fn = gen_fname(tag, seed)
        fn = "test_QSA_$(tag)_tau$(τ)_seed$(seed).2.txt"
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
            accrate[1] = λ * accrate[1] + (1 - λ) * acc / it
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    userrr = false

    fst = true; t0 = 0.0; C = nothing; it0 = 0
    C = RRRMC.Config(N * M)
    for k = 1:M
        copy!(C.s, (k-1)*N + 1, C1.s, 1, N)
    end
    local X
    for Γ = linspace(Γ0, Γ1, num_steps)
    #Γ = Γ0
    #for s = 1:num_steps
        info("Γ=$Γ β=$β ar=$(accrate[1]) useRRR=$userrr")
        sMC = (tag == :rrr && userrr) ? RRRMC.rrrMC : RRRMC.standardMC
        hook, cleanup, t0 = gen_hook(string(tag), β, Γ, fst, t0, it0)
        try
            X = RRRMC.GraphQPercT(Xs, M, Γ, β)
            @time E, C = sMC(X, β, iters, seed=sseed, step=itst, hook=hook, C0=C)
        finally
            cleanup()
        end
        fst = false
        it0 += iters
        #sseed += 624234
        sseed = 0
        accrate[1] < accrthresh && (userrr = true)

        β += (β1 - β) * 0.1
        #Γ += (Γ1 - Γ) * 0.1
    end

    Es = RRRMC.Renergies(X)
    @show Es
    @show minimum(Es)

    return C, Xs, X
end

findb(b, it, ns, x) = abs(b * (-1 + (b / (b + it))^(-1/ns)) - x)

function test_QSA_alt(;N::Integer = 1001, α::Real = 0.6, M::Integer = 32, alg = :rrr, τ::Integer = 10^2, seed = 78821000027346, num_steps::Integer = 30, tag = "")
    β = 32.0
    β1 = 32.0

    C1, Xs = test_SA(T0 = 3.0, T1 = 1/β, τ = 100, N = N, α = α, seed = seed + 78123472837432554)
    N = RRRMC.getN(Xs)
    P = Xs.P

    Γ0 = 10.0
    Γ1 = 1e-2

    tot_iters = τ * N

    #num_steps = 30

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

    @assert τ ≥ 100
    λ = (1 - 100/τ)^10 # very empirical formula...
    accrthresh = 1.0 # 0.07

    sseed = seed

    accrate = Float64[0.1]
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
            accrate[1] = λ * accrate[1] + (1 - λ) * acc / it
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, t0
    end

    userrr = false

    fst = true; t0 = 0.0; C = nothing; it0 = 0
    C = RRRMC.Config(N * M)
    for k = 1:M
        copy!(C.s, (k-1)*N + 1, C1.s, 1, N)
    end
    local X
    for s = 1:num_steps
        info("Γ=$(Γ[s]) β=$β ar=$(accrate[1]) useRRR=$userrr")
        sMC = (alg == :rrr && userrr) ? RRRMC.rrrMC : RRRMC.standardMC
        hook, cleanup, t0 = gen_hook(string(alg), tag, β, Γ[s], fst, t0, it0)
        try
            X = RRRMC.GraphQPercT(Xs, M, Γ[s], β)
            @time E, C = sMC(X, β, iters[s], seed=sseed, step=itst[s], hook=hook, C0=C)
        finally
            cleanup()
        end
        fst = false
        it0 += iters[s]
        #sseed += 624234
        sseed = 0
        accrate[1] < accrthresh && (userrr = true)

        β += (β1 - β) * 0.1
    end

    Es = RRRMC.Renergies(X)
    @show Es
    @show minimum(Es)

    return C, Xs, X
end


end # module
