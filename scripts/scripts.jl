## The scripts used to produce the results in the paper
## "A method to reduce the rejection rate in Monte Carlo Markov Chains
## on Ising spin models"

module Scripts

using JLD
using Compat

include("../src/RRRMC.jl")

function to_mat(Cv::Vector{BitVector})
    samples = length(Cv)
    N = length(Cv[1])
    Cs = BitMatrix(N, samples)
    for i = 1:samples
        Cs[:,i] = Cv[i]
    end
    return Cs
end

function test_RRG(; N = 10_000,
                    K = 3,
                    β = 2.0,
                    iters = 10^14,
                    step = 10^4,
                    seedx = 8370000274,
                    seed0 = 6540000789,
                    seedst = 10_000,
                    seedstx = nothing,
                    ntests = 10,
                    force = false,
                    met_factor = 3.7,  # β=3 → 4.0    β=4 → 4.5
                    bkl_factor = 94.9, # β=3 → 768.6  β=4 → 6082.7
                    rrr_factor = 1.0,
                    wtm_factor = 53.0, # β=3 → 412.1  β=4 → 3375.2
                    t_limit = 40.0,
                    tag = "",
                    algs = [:met, :bkl, :rrr, :wtm]
                 )

    @assert all(a ∈ [:met, :bkl, :rrr, :wtm] for a in algs)

    seedstx::Int = (seedstx ≡ nothing ? seedst : seedstx)

    dirname = "output_RRG_N$(N)_K$(K)_beta$(β)_tmax$(t_limit)_step$(step)$(tag)"
    gen_fname(alg, seed) = joinpath(dirname, "output_$(alg)_sx$(seedx)_s$(seed).txt")
    gen_Cfname(alg, seedx, seed) = joinpath(dirname, "Cs_$(alg)_sx$(seedx)_s$(seed).jld")

    samples = iters ÷ step

    function gen_hook(alg, seed)
        isdir(dirname) || mkdir(dirname)
        fn = gen_fname(alg, seed)
        !force && isfile(fn) && error("file $fn exists")
        f = open(fn, "w")
        Cv = Vector{BitVector}() #BitArray(N, samples)
        t0 = time()
        println(f, "#mctime acc E clocktime")
        hook = (mct, X, C, acc, E) -> begin
            t = time() - t0
            push!(Cv, copy(C.s))
            println(f, "$mct $acc $E $t")
            return t < t_limit
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, Cv
    end

    srand(seedx)
    X = RRRMC.GraphRRG(N, K)

    info("compile...")
    RRRMC.standardMC(X, β, 10^3)
    RRRMC.bklMC(X, β, 10^3)
    RRRMC.wtmMC(X, β, 10^3)
    RRRMC.rrrMC(X, β, 10^3)

    seed = seed0
    for tst = 1:ntests
        info("### SEED = $seed SEEDx = $seedx")

        srand(seedx)
        X = RRRMC.GraphRRG(N, K)

        if :met in algs
            info("# Metropolis")
            gc()
            rstep = round(Int, step * met_factor)
            riters = rstep * samples
            hook, cleanup, met_Cv = gen_hook("met", seed)
            try
                @time met_Es, met_C = RRRMC.standardMC(X, β, riters, step=rstep, seed=seed, hook=hook)
            finally
                cleanup()
            end
            met_Cs = to_mat(met_Cv)
            save(gen_Cfname("met", seedx, seed), Dict("Cs"=>met_Cs))
        end
        if :bkl in algs
            info("# BKL")
            gc()
            rstep = round(Int, step * bkl_factor)
            riters = rstep * samples
            hook, cleanup, bkl_Cv = gen_hook("bkl", seed)
            try
                @time bkl_Es, bkl_C = RRRMC.bklMC(X, β, riters, step=rstep, seed=seed, hook=hook)
            finally
                cleanup()
            end
            bkl_Cs = to_mat(bkl_Cv)
            save(gen_Cfname("bkl", seedx, seed), Dict("Cs"=>bkl_Cs))
        end
        if :rrr in algs
            info("# RRR")
            gc()
            rstep = round(Int, step * rrr_factor)
            riters = rstep * samples
            hook, cleanup, rrr_Cv = gen_hook("rrr", seed)
            try
                @time rrr_Es, rrr_C = RRRMC.rrrMC(X, β, riters, step=rstep, seed=seed, hook=hook)
            finally
                cleanup()
            end
            rrr_Cs = to_mat(rrr_Cv)
            save(gen_Cfname("rrr", seedx, seed), Dict("Cs"=>rrr_Cs))
        end
        if :wtm in algs
            info("# WTM")
            gc()
            rtstep = step * wtm_factor
            hook, cleanup, wtm_Cv = gen_hook("wtm", seed)
            try
                @time wtm_Es, wtm_C = RRRMC.wtmMC(X, β, samples, step=rtstep, seed=seed, hook=hook)
            finally
                cleanup()
            end
            wtm_Cs = to_mat(wtm_Cv)
            save(gen_Cfname("wtm", seedx, seed), Dict("Cs"=>wtm_Cs))
        end

        seed += seedst
        seedx += seedstx

        println(STDERR)
    end
end

function test_RRGCont(; N = 10_000,
                        K = 3,
                        β = 2.0,
                        iters = 10^14,
                        step = 10^4,
                        seedx = 8370000274,
                        seed0 = 6540000789,
                        seedst = 10_000,
                        seedstx = nothing,
                        ntests = 10,
                        force = false,
                        met_factor = 8.0,  # β=3 → 7.3  β=4 → 7.5
                        bkl_factor = 16.5, # β=3 → 32.8 β=4 → 46.3
                        rrr_factor = 1.0,
                        wtm_factor = 20.5, # β=3 → 38.0 β=4 → 57.2
                        t_limit = 40.0,
                        tag = "",
                        algs = [:met, :bkl, :rrr, :wtm]
                     )

    @assert all(a ∈ [:met, :bkl, :rrr, :wtm] for a in algs)

    seedstx::Int = (seedstx ≡ nothing ? seedst : seedstx)

    dirname = "output_RRGCont_N$(N)_K$(K)_beta$(β)_tmax$(t_limit)_step$(step)$(tag)"
    gen_fname(alg, seed) = joinpath(dirname, "output_$(alg)_sx$(seedx)_s$(seed).txt")
    gen_Cfname(alg, seedx, seed) = joinpath(dirname, "Cs_$(alg)_sx$(seedx)_s$(seed).jld")

    samples = iters ÷ step

    function gen_hook(alg, seed)
        isdir(dirname) || mkdir(dirname)
        fn = gen_fname(alg, seed)
        !force && isfile(fn) && error("file $fn exists")
        f = open(fn, "w")
        Cv = Vector{BitVector}() #BitArray(N, samples)
        t0 = time()
        println(f, "#mctime acc E clocktime")
        hook = (mct, X, C, acc, E) -> begin
            t = time() - t0
            push!(Cv, copy(C.s))
            println(f, "$mct $acc $E $t")
            return t < t_limit
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, Cv
    end

    srand(seedx)
    X = RRRMC.GraphRRGNormal(N, K)

    info("compile...")
    RRRMC.standardMC(X, β, 10^3)
    RRRMC.bklMC(X, β, 10^3)
    RRRMC.wtmMC(X, β, 10^3)
    RRRMC.rrrMC(X, β, 10^3)

    seed = seed0
    for tst = 1:ntests
        info("### SEED = $seed SEEDx = $seedx")

        srand(seedx)
        X = RRRMC.GraphRRGNormal(N, K)

        if :met in algs
            info("# Metropolis")
            gc()
            rstep = round(Int, step * met_factor)
            riters = rstep * samples
            hook, cleanup, met_Cv = gen_hook("met", seed)
            try
                @time met_Es, met_C = RRRMC.standardMC(X, β, riters, step=rstep, seed=seed, hook=hook)
            finally
                cleanup()
            end
            met_Cs = to_mat(met_Cv)
            save(gen_Cfname("met", seedx, seed), Dict("Cs"=>met_Cs))
        end
        if :bkl in algs
            info("# BKL")
            gc()
            rstep = round(Int, step * bkl_factor)
            riters = rstep * samples
            hook, cleanup, bkl_Cv = gen_hook("bkl", seed)
            try
                @time bkl_Es, bkl_C = RRRMC.bklMC(X, β, riters, step=rstep, seed=seed, hook=hook)
            finally
                cleanup()
            end
            bkl_Cs = to_mat(bkl_Cv)
            save(gen_Cfname("bkl", seedx, seed), Dict("Cs"=>bkl_Cs))
        end
        if :rrr in algs
            info("# RRR")
            gc()
            rstep = round(Int, step * rrr_factor)
            riters = rstep * samples
            hook, cleanup, rrr_Cv = gen_hook("rrr", seed)
            try
                @time rrr_Es, rrr_C = RRRMC.rrrMC(X, β, riters, step=rstep, seed=seed, hook=hook)
            finally
                cleanup()
            end
            rrr_Cs = to_mat(rrr_Cv)
            save(gen_Cfname("rrr", seedx, seed), Dict("Cs"=>rrr_Cs))
        end
        if :wtm in algs
            info("# WTM")
            gc()
            rtstep = step * wtm_factor
            hook, cleanup, wtm_Cv = gen_hook("wtm", seed)
            try
                @time wtm_Es, wtm_C = RRRMC.wtmMC(X, β, samples, step=rtstep, seed=seed, hook=hook)
            finally
                cleanup()
            end
            wtm_Cs = to_mat(wtm_Cv)
            save(gen_Cfname("wtm", seedx, seed), Dict("Cs"=>wtm_Cs))
        end

        seed += seedst
        seedx += seedstx

        println(STDERR)
    end
end

# same as (2a-1) ⋅ (2b-1)
# note: no length checks
function pm1dot(a::BitVector, b::BitVector)
    # one way to write it avoiding allocations:
    # 4 * (a ⋅ b) - 2sum(a) - 2sum(b) + length(a)

    # ugly but slightly faster (length(a)-2sum(a ⊻ b) without allocations):
    l = length(a)
    ac = a.chunks
    bc = b.chunks
    @inbounds @simd for i = 1:length(ac)
        l -= 2 * count_ones(ac[i] ⊻ bc[i])
    end
    return l
end

function parsets(outfname, tag::Symbol)
    info("parse $outfname")
    # ratio = 1.0
    # if tag == :wtm
    #     final_t = 0.0
    #     final_gt = 1.0
    #     ratio = open(outfname) do f
    #         for l in eachline(f)
    #             eof(f) || continue
    #             sl = split(l)
    #             final_t = parse(Float64, sl[4])
    #             final_gt = parse(Float64, sl[1])
    #         end
    #         return final_t / final_gt
    #     end
    # end

    ts = open(outfname) do of
        l = readline(of)
        @assert startswith(l, "#")
        ts = Float64[]
        for l in eachline(of)
            @assert !startswith(l, "#")
            sl = split(l)
            @assert length(sl) ≥ 4
            #t = tag == :wtm ? parse(Float64, sl[1]) : parse(Float64, sl[4])
            t = parse(Float64, sl[4])
            push!(ts, t)
        end
        return ts
    end

    return ts #, ratio
end

type LogRange{T}
    i0::T
    i1::T
    st0::Float64
    incr::Float64
    function LogRange(i0::T, i1::T, st0::Real, incr::Real)
        @assert st0 > 0
        @assert incr ≥ 1
        return new(i0, i1, st0, incr)
    end
end

LogRange{T}(i0::T, i1::T, st0::Real, incr::Real) = LogRange{T}(i0, i1, st0, incr)

const LRINCR = 1.5

LogRange(l::Int) = LogRange(1, l, 1, LRINCR)

Base.start(lr::LogRange) = (lr.i0, lr.st0)
Base.done{T}(lr::LogRange, stat::Tuple{T,Float64}) = stat[1] > lr.i1
function Base.next{T<:Integer}(lr::LogRange{T}, stat::Tuple{T,Float64})
    i, st = stat
    inext = i + max(1, round(T, st))
    stnext = st * lr.incr
    return i, (inext, stnext)
end
function Base.next{T}(lr::LogRange{T}, stat::Tuple{T,Float64})
    i, st = stat
    inext = i + st
    stnext = st * lr.incr
    return i, (inext, stnext)
end

trunclen!(v::Vector, i::Integer) = deleteat!(v, (i+1):length(v))
trunclen!{N}(vs::NTuple{N,Vector}, i::Integer) = for v in vs; trunclen!(v, i); end
trunclen!(vs::Vector...) = trunclen!(vs, minimum(length(v) for v in vs))

function get_ts_range(tsm::Vector{Float64}, t0::Real)
    i = findfirst(t->t ≥ t0, tsm)
    j = findfirst(t->t ≥ 2t0, tsm)
    return i, j
end

function parseovs(Cfname::AbstractString, tsm::Vector{Float64}, lr::LogRange{Float64})
    info("parse $Cfname")
    Cs::BitMatrix = load(Cfname, "Cs")
    N, l = size(Cs)

    Cv = BitVector[Cs[:,i] for i = 1:l]

    mq2s = Float64[]
    sq2s = Float64[]

    for t_st in lr
        #println("   i=$i/$l")
        i, j = get_ts_range(tsm, t_st)
        i == 0 && break
        j == 0 && (j = length(tsm)+1)

        mq2 = 0.0
        mq4 = 0.0
        n = 0

        for i1 = i:(j-2), j1 = (i+1):(j-1)
            si = Cv[i1]
            sj = Cv[j1]
            q2 = (pm1dot(si, sj) / N)^2
            mq2 += q2
            mq4 += q2^2
            n += 1
        end

        mq2 /= n
        mq4 /= n
        sq2 = √max(0.0, mq4 - mq2^2) # TODO: Bessel's correction?

        push!(mq2s, mq2)
        push!(sq2s, sq2)
    end
    return mq2s, sq2s
end

function parsexovs(Cfname1::AbstractString, Cfname2::AbstractString, ts1::Vector{Float64}, ts2::Vector{Float64}, lr::LogRange{Float64})
    info("parse $Cfname1 + $Cfname2")
    Cs1::BitMatrix = load(Cfname1, "Cs")
    Cs2::BitMatrix = load(Cfname2, "Cs")
    N, l1 = size(Cs1)
    l2 = size(Cs2, 2)
    @assert size(Cs2, 1) == N
    #l = min(l1, l2)

    Cv1 = BitVector[Cs1[:,i] for i = 1:l1]
    Cv2 = BitVector[Cs2[:,i] for i = 1:l2]

    mx2s = Float64[]
    sx2s = Float64[]

    for t_st in lr

        i1, j1 = get_ts_range(ts1, t_st)
        i2, j2 = get_ts_range(ts2, t_st)

        (i1 == 0 || i2 == 0) && break

        j1 == 0 && (j1 = length(ts1)+1)
        j2 == 0 && (j2 = length(ts2)+1)

        mx2 = 0.0
        mx4 = 0.0
        n = 0

        for k1 = i1:(j1-1), k2 = i2:(j2-1)

            s1 = Cv1[k1]
            s2 = Cv2[k2]

            x2 = (pm1dot(s1, s2) / N)^2
            #push!(x2s, x2)
            mx2 += x2
            mx4 += x2^2
            n += 1
        end

        mx2 /= n
        mx4 /= n

        sx2 = √max(0.0, mx4 - mx2^2) # TODO: Bessel's correction?

        push!(mx2s, mx2)
        push!(sx2s, sx2)
    end
    return mx2s, sx2s
end

function stats_overlaps(dirname::AbstractString, tags::Vector{Symbol}, sx::Int = 8370000274, s1::Int = 6540000789, s2::Int = 5430000678, step::Float64 = 30.0/2^8, incr::Float64 = 2.0)
    for tag in tags
        stats_overlaps(dirname, tag, sx, s1, s2, step, incr)
    end
end

function stats_overlaps(dirname::AbstractString, tag::Symbol, sx::Int, s1::Int, s2::Int, step::Float64 = 30.0/2^8, incr::Float64 = 2.0)
    isdir(dirname) || error("directory $dirname not found")
    tags = [:met, :bkl, :rrr, :wtm]
    @assert tag ∈ tags

    outfname1 = joinpath(dirname, "output_$(tag)_sx$(sx)_s$(s1).txt")
    outfname2 = joinpath(dirname, "output_$(tag)_sx$(sx)_s$(s2).txt")
    Cfname1 = joinpath(dirname, "Cs_$(tag)_sx$(sx)_s$(s1).jld")
    Cfname2 = joinpath(dirname, "Cs_$(tag)_sx$(sx)_s$(s2).jld")

    isfile(outfname1) || error("file $outfname1 not found")
    isfile(outfname2) || error("file $outfname2 not found")
    isfile(Cfname1) || error("file $Cfname1 not found")
    isfile(Cfname2) || error("file $Cfname2 not found")

    fname = joinpath(dirname, "overlaps_$(tag)_sx$(sx).txt")

    lr = LogRange(step, 10000.0, step, incr)

    ts1 = parsets(outfname1, tag)
    ts2 = parsets(outfname2, tag)

    #l = min(length(ts1), length(ts2))

    #trunclen!(ts1, ts2)

    mq2s1, sq2s1 = parseovs(Cfname1, ts1, lr)
    mq2s2, sq2s2 = parseovs(Cfname2, ts2, lr)

    trunclen!(mq2s1, sq2s1, mq2s2, sq2s2)

    mq2s = (mq2s1 + mq2s2) / 2
    sq2s = (sq2s1 + sq2s2) / 2

    mx2s, sx2s = parsexovs(Cfname1, Cfname2, ts1, ts2, lr)

    r = min(length(mx2s), length(mq2s))

    info("write output")
    open(fname, "w") do f
        println(f, "#time self std cross std")
        for (i,t) = enumerate(lr)
            i > r && break
            mq2, sq2 = mq2s[i], sq2s[i]
            mx2, sx2 = mx2s[i], sx2s[i]

            (isnan(mq2) || isnan(sq2) || isnan(mx2) || isnan(sx2)) && (warn("Found NaN at t=$t (mq2=$mq2, sq2=$sq2, mx2=$mx2, sx2=$sx2)"); continue)

            println(f, "$t $mq2 $sq2 $mx2 $sx2")
        end
    end
    info("done")

    return
end

function stats_overlaps_all(dirname::AbstractString; outlfrac::Float64 = 0.0)
    isdir(dirname) || error("directory $dirname not found")
    tags = [:met, :bkl, :rrr, :wtm]

    diffs = Dict{Symbol, Dict{AbstractString,Float64}}()
    files = readdir(dirname)
    for fn in files
        isfile(joinpath(dirname, fn)) || continue
        startswith(fn, "stats_") && continue
        startswith(fn, "output_") && continue
        startswith(fn, "Cs_") && continue
        ismatch(r"^overlaps_(...)_sx\d+\.txt$", fn) || continue
        tag = Symbol(match(r"overlaps_(...)_sx\d+\.txt", fn).captures[1])
        @assert tag ∈ tags
        dt = get!(diffs, tag, Dict{AbstractString,Float64}())
        @assert !haskey(dt, fn)
        open(joinpath(dirname, fn)) do f
            for l in eachline(f)
                eof(f) || continue
                sl = split(l)
                @assert length(sl) ≥ 4
                q2m = parse(Float64, sl[2])
                x2m = parse(Float64, sl[4])
                dt[fn] = q2m - x2m
            end
        end
    end
    outliers = Dict{Symbol, Set{AbstractString}}()
    for (tag, dt) in diffs
        sd = sort!([(fn,d) for (fn,d) in dt], by=x->x[2], rev=true)
        sd = sd[1:round(Int, length(sd)*outlfrac)]
        outliers[tag] = Set{AbstractString}([s[1] for s in sd])
    end

    fname = joinpath(dirname, "stats_overlaps.txt")
    open(fname, "w") do sf
        println(sf, "#time ", join(("$(tag)_self std $(tag)_cross std" for tag in tags), " "))
        q2mdict = Dict{Symbol, Dict{Float64, Vector{Float64}}}(t => Dict{Float64, Vector{Float64}}() for t in tags)
        q2sdict = Dict{Symbol, Dict{Float64, Vector{Float64}}}(t => Dict{Float64, Vector{Float64}}() for t in tags)
        x2mdict = Dict{Symbol, Dict{Float64, Vector{Float64}}}(t => Dict{Float64, Vector{Float64}}() for t in tags)
        x2sdict = Dict{Symbol, Dict{Float64, Vector{Float64}}}(t => Dict{Float64, Vector{Float64}}() for t in tags)

        files = readdir(dirname)
        for fn in files
            isfile(joinpath(dirname, fn)) || continue
            startswith(fn, "stats_") && continue
            startswith(fn, "output_") && continue
            startswith(fn, "Cs_") && continue
            ismatch(r"^overlaps_(...)_sx\d+\.txt$", fn) || (warn("skipping file $fn"); continue)
            tag = Symbol(match(r"overlaps_(...)_sx\d+\.txt", fn).captures[1])
            @assert tag ∈ tags
            fn ∈ outliers[tag] && (info("skip $fn"); continue)
            q2md_tag = q2mdict[tag]
            q2sd_tag = q2sdict[tag]
            x2md_tag = x2mdict[tag]
            x2sd_tag = x2sdict[tag]
            open(joinpath(dirname, fn)) do f
                for l in eachline(f)
                    ismatch(r"^\s*#", l) && continue
                    sl = split(l)
                    @assert length(sl) ≥ 4
                    t_st = parse(Float64, sl[1])
                    q2m = parse(Float64, sl[2])
                    q2s = parse(Float64, sl[3])
                    x2m = parse(Float64, sl[4])
                    x2s = parse(Float64, sl[5])

                    q2md_vec = get!(q2md_tag, t_st, Float64[])
                    q2sd_vec = get!(q2sd_tag, t_st, Float64[])
                    x2md_vec = get!(x2md_tag, t_st, Float64[])
                    x2sd_vec = get!(x2sd_tag, t_st, Float64[])
                    push!(q2md_vec, q2m)
                    push!(q2sd_vec, q2s)
                    push!(x2md_vec, x2m)
                    push!(x2sd_vec, x2s)
                end
            end
        end

        all_t_st = sort!(union([collect(keys(q2md_tag)) for q2md_tag in values(q2mdict)]...))
        L = length(all_t_st)

        hasnans = false

        mq2mdict = Dict{Symbol, Vector{Float64}}(tag => zeros(L) for tag in tags)
        mq2sdict = Dict{Symbol, Vector{Float64}}(tag => zeros(L) for tag in tags)
        mx2mdict = Dict{Symbol, Vector{Float64}}(tag => zeros(L) for tag in tags)
        mx2sdict = Dict{Symbol, Vector{Float64}}(tag => zeros(L) for tag in tags)

        for tag in tags
            q2md_tag = q2mdict[tag]
            x2md_tag = x2mdict[tag]

            mq2md_tag = mq2mdict[tag]
            mq2sd_tag = mq2sdict[tag]
            mx2md_tag = mx2mdict[tag]
            mx2sd_tag = mx2sdict[tag]

            for i = 1:L
                t_st = all_t_st[i]
                for (srcd, md, sd) in [(q2md_tag, mq2md_tag, mq2sd_tag),
                                       (x2md_tag, mx2md_tag, mx2sd_tag)]
                    if !haskey(srcd, t_st)
                        hasnans = true
                        md[i] = NaN
                        sd[i] = NaN
                        continue
                    end
                    md[i] = mean(srcd[t_st])
                    sd[i] = std(srcd[t_st]) / √length(srcd[t_st])
                end
            end
        end

        hasnans && warn("empty bins found (increase step?)")

        for i = 1:L
            t_st = all_t_st[i]
            println(sf, t_st, " ", join(vec([string(md[tag][i]) for md in [mq2mdict, mq2sdict, mx2mdict, mx2sdict], tag in tags]), " "))
        end
    end
end

function stats_overlaps_all_diff(dirname::AbstractString; outlfrac::Float64 = 0.0)
    isdir(dirname) || error("directory $dirname not found")
    tags = [:met, :bkl, :rrr, :wtm]

    diffs = Dict{Symbol, Dict{AbstractString,Float64}}()
    files = readdir(dirname)
    for fn in files
        isfile(joinpath(dirname, fn)) || continue
        startswith(fn, "stats_") && continue
        startswith(fn, "output_") && continue
        startswith(fn, "Cs_") && continue
        ismatch(r"^overlaps_(...)_sx\d+\.txt$", fn) || (warn("skipping file $fn"); continue)
        tag = Symbol(match(r"overlaps_(...)_sx\d+\.txt", fn).captures[1])
        @assert tag ∈ tags
        dt = get!(diffs, tag, Dict{AbstractString,Float64}())
        @assert !haskey(dt, fn)
        open(joinpath(dirname, fn)) do f
            for l in eachline(f)
                eof(f) || continue
                sl = split(l)
                @assert length(sl) ≥ 4
                q2m = parse(Float64, sl[2])
                x2m = parse(Float64, sl[4])
                dt[fn] = q2m - x2m
            end
        end
    end
    outliers = Dict{Symbol, Set{AbstractString}}()
    for (tag, dt) in diffs
        sd = sort!([(fn,d) for (fn,d) in dt], by=x->x[2], rev=true)
        sd = sd[1:round(Int, length(sd)*outlfrac)]
        outliers[tag] = Set{AbstractString}([s[1] for s in sd])
    end
    @show outliers

    fname = joinpath(dirname, "stats_overlaps_diff.txt")
    open(fname, "w") do sf
        println(sf, "#time ", join(("$(tag)_diff std" for tag in tags), " "))
        q2mdict = Dict{Symbol, Dict{Float64, Vector{Float64}}}(t => Dict{Float64, Vector{Float64}}() for t in tags)
        q2sdict = Dict{Symbol, Dict{Float64, Vector{Float64}}}(t => Dict{Float64, Vector{Float64}}() for t in tags)
        x2mdict = Dict{Symbol, Dict{Float64, Vector{Float64}}}(t => Dict{Float64, Vector{Float64}}() for t in tags)
        x2sdict = Dict{Symbol, Dict{Float64, Vector{Float64}}}(t => Dict{Float64, Vector{Float64}}() for t in tags)

        files = readdir(dirname)
        for fn in files
            isfile(joinpath(dirname, fn)) || continue
            startswith(fn, "stats_") && continue
            startswith(fn, "output_") && continue
            startswith(fn, "Cs_") && continue
            ismatch(r"^overlaps_(...)_sx\d+\.txt$", fn) || (warn("skipping file $fn"); continue)
            tag = Symbol(match(r"overlaps_(...)_sx\d+\.txt", fn).captures[1])
            @assert tag ∈ tags
            fn ∈ outliers[tag] && (info("skip $fn"); continue)
            q2md_tag = q2mdict[tag]
            q2sd_tag = q2sdict[tag]
            x2md_tag = x2mdict[tag]
            x2sd_tag = x2sdict[tag]
            open(joinpath(dirname, fn)) do f
                for l in eachline(f)
                    ismatch(r"^\s*#", l) && continue
                    sl = split(l)
                    @assert length(sl) ≥ 4
                    t_st = parse(Float64, sl[1])
                    q2m = parse(Float64, sl[2])
                    q2s = parse(Float64, sl[3])
                    x2m = parse(Float64, sl[4])
                    x2s = parse(Float64, sl[5])

                    q2md_vec = get!(q2md_tag, t_st, Float64[])
                    q2sd_vec = get!(q2sd_tag, t_st, Float64[])
                    x2md_vec = get!(x2md_tag, t_st, Float64[])
                    x2sd_vec = get!(x2sd_tag, t_st, Float64[])
                    push!(q2md_vec, q2m)
                    push!(q2sd_vec, q2s)
                    push!(x2md_vec, x2m)
                    push!(x2sd_vec, x2s)
                end
            end
        end

        all_t_st = sort!(union([collect(keys(q2md_tag)) for q2md_tag in values(q2mdict)]...))
        L = length(all_t_st)

        hasnans = false

        md2mdict = Dict{Symbol, Vector{Float64}}(tag => zeros(L) for tag in tags)
        md2sdict = Dict{Symbol, Vector{Float64}}(tag => zeros(L) for tag in tags)

        for tag in tags
            q2md_tag = q2mdict[tag]
            x2md_tag = x2mdict[tag]

            md2md_tag = md2mdict[tag]
            md2sd_tag = md2sdict[tag]

            for i = 1:L
                t_st = all_t_st[i]
                if !haskey(q2md_tag, t_st) || !haskey(x2md_tag, t_st)
                    hasnans = true
                    md2md_tag[i] = NaN
                    md2sd_tag[i] = NaN
                    continue
                end
                dv = q2md_tag[t_st] - x2md_tag[t_st]
                #@show tag, t_st, sum(dv .< 0)
                md2md_tag[i] = mean(dv)
                md2sd_tag[i] = std(dv) / √length(dv)
            end
        end

        hasnans && warn("empty bins found (increase step?)")

        for i = 1:L
            t_st = all_t_st[i]
            println(sf, t_st, " ", join(vec([string(md[tag][i]) for md in [md2mdict, md2sdict], tag in tags]), " "))
        end
    end
end

function test_QIsing(; N = 1_024,
                       M = 16,
                       β = 2.0,
                       Γ = 0.3,
                       iters = 10^14,
                       step = 10^4,
                       seedx = 8370000274,
                       seed0 = 6540000789,
                       seedst = 10_000,
                       seedstx = nothing,
                       ntests = 10,
                       force = false,
                       met_factor = 15.74,
                       rrr_factor = 1.0,
                       t_limit = 250.0,
                       tag = "",
                       algs = [:met, :rrr]
                    )

    @assert all(a ∈ [:met, :rrr] for a in algs)

    seedstx::Int = (seedstx ≡ nothing ? seedst : seedstx)

    dirname = "output_QIsing_N$(N)_M$(M)_beta$(β)_Gamma$(Γ)_tmax$(t_limit)_step$(step)$(tag)"
    gen_fname(alg, seed) = joinpath(dirname, "output_$(alg)_sx$(seedx)_s$(seed).txt")
    gen_Cfname(alg, seedx, seed) = joinpath(dirname, "Cs_$(alg)_sx$(seedx)_s$(seed).jld")

    samples = iters ÷ step

    function gen_hook(alg, seed)
        isdir(dirname) || mkdir(dirname)
        fn = gen_fname(alg, seed)
        !force && isfile(fn) && error("file $fn exists")
        f = open(fn, "w")
        Cv = Vector{BitVector}() #BitArray(N, samples)
        t0 = time()
        println(f, "#mctime acc QE clocktime")
        hook = (mct, X, C, acc, E) -> begin
            t = time() - t0
            push!(Cv, copy(C.s))
            QE = RRRMC.Qenergy(X, C)
            println(f, "$mct $acc $QE $t")
            return t < t_limit
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, Cv
    end

    srand(seedx)
    X = RRRMC.GraphQSKT(N, M, Γ, β)

    info("compile...")
    RRRMC.standardMC(X, β, 10^3)
    RRRMC.rrrMC(X, β, 10^3)

    seed = seed0
    for tst = 1:ntests
        info("### SEED = $seed SEEDx = $seedx")

        srand(seedx)
        X = RRRMC.GraphQSKT(N, M, Γ, β)

        if :met in algs
            info("# Metropolis")
            gc()
            rstep = round(Int, step * met_factor)
            riters = rstep * samples
            hook, cleanup, met_Cv = gen_hook("met", seed)
            try
                @time met_Es, met_C = RRRMC.standardMC(X, β, riters, step=rstep, seed=seed, hook=hook)
            finally
                cleanup()
            end
            met_Cs = to_mat(met_Cv)
            save(gen_Cfname("met", seedx, seed), Dict("Cs"=>met_Cs))
        end
        if :rrr in algs
            info("# RRR")
            gc()
            rstep = round(Int, step * rrr_factor)
            riters = rstep * samples
            hook, cleanup, rrr_Cv = gen_hook("rrr", seed)
            try
                @time rrr_Es, rrr_C = RRRMC.rrrMC(X, β, riters, step=rstep, seed=seed, hook=hook)
            finally
                cleanup()
            end
            rrr_Cs = to_mat(rrr_Cv)
            save(gen_Cfname("rrr", seedx, seed), Dict("Cs"=>rrr_Cs))
        end

        seed += seedst
        seedx += seedstx

        println(STDERR)
    end
end

function test_REIsing(; N = 1_024,
                        M = 5,
                        β = 0.4,
                        γ = 2.0,
                        iters = 10^14,
                        step = 10^4,
                        seedx = 8370000274,
                        seed0 = 6540000789,
                        seedst = 10_000,
                        seedstx = nothing,
                        ntests = 10,
                        force = false,
                        met_factor = 20.8, # γ=3 → 24.6 # γ=4 → 13.9 # γ=5 → 6.4
                        rrr_factor = 1.0,
                        t_limit = 250.0,
                        tag = "",
                        algs = [:met, :rrr]
                     )

    @assert all(a ∈ [:met, :rrr] for a in algs)

    seedstx::Int = (seedstx ≡ nothing ? seedst : seedstx)

    dirname = "output_REIsing_N$(N)_M$(M)_beta$(β)_gamma$(γ)_tmax$(t_limit)_step$(step)$(tag)"
    gen_fname(alg, seed) = joinpath(dirname, "output_$(alg)_sx$(seedx)_s$(seed).txt")
    gen_Cfname(alg, seedx, seed) = joinpath(dirname, "Cs_$(alg)_sx$(seedx)_s$(seed).jld")

    samples = iters ÷ step

    function gen_hook(alg, seed)
        isdir(dirname) || mkdir(dirname)
        fn = gen_fname(alg, seed)
        !force && isfile(fn) && error("file $fn exists")
        f = open(fn, "w")
        Cv = Vector{BitVector}() #BitArray(N, samples)
        t0 = time()
        println(f, "#mctime acc meanRE clocktime E")
        hook = (mct, X, C, acc, E) -> begin
            t = time() - t0
            push!(Cv, copy(C.s))
            meanRE = mean(RRRMC.REenergies(X))
            println(f, "$mct $acc $meanRE $t $E")
            return t < t_limit
        end
        cleanup = () -> begin
            close(f)
        end
        return hook, cleanup, Cv
    end

    srand(seedx)
    X = RRRMC.GraphSKRE(N, M, γ, β)

    info("compile...")
    RRRMC.standardMC(X, β, 10^3)
    RRRMC.rrrMC(X, β, 10^3)

    seed = seed0
    for tst = 1:ntests
        info("### SEED = $seed SEEDx = $seedx")

        srand(seedx)
        X = RRRMC.GraphSKRE(N, M, γ, β)

        if :met in algs
            info("# Metropolis")
            gc()
            rstep = round(Int, step * met_factor)
            riters = rstep * samples
            hook, cleanup, met_Cv = gen_hook("met", seed)
            try
                @time met_Es, met_C = RRRMC.standardMC(X, β, riters, step=rstep, seed=seed, hook=hook)
            finally
                cleanup()
            end
            met_Cs = to_mat(met_Cv)
            save(gen_Cfname("met", seedx, seed), Dict("Cs"=>met_Cs))
        end
        if :rrr in algs
            info("# RRR")
            gc()
            rstep = round(Int, step * rrr_factor)
            riters = rstep * samples
            hook, cleanup, rrr_Cv = gen_hook("rrr", seed)
            try
                @time rrr_Es, rrr_C = RRRMC.rrrMC(X, β, riters, step=rstep, seed=seed, hook=hook)
            finally
                cleanup()
            end
            rrr_Cs = to_mat(rrr_Cv)
            save(gen_Cfname("rrr", seedx, seed), Dict("Cs"=>rrr_Cs))
        end

        seed += seedst
        seedx += seedstx

        println(STDERR)
    end
end
function stats_time(dirname::AbstractString; step::Float64 = 2.0, algs = [:met, :bkl, :rrr, :wtm])
    isdir(dirname) || error("directory $dirname not found")

    @assert all(a ∈ [:met, :bkl, :rrr, :wtm] for a in algs)

    fname = joinpath(dirname, "stats_time.txt")

    #tname = joinpath(dirname, "test.tmp.txt")
    #tf = open(tname, "w")

    open(fname, "w") do sf
        println(sf, "#time ", join((string(tag, " std") for tag in algs), " "))
        Edict = Dict{Symbol, Dict{Float64, Vector{Float64}}}(t => Dict{Float64, Vector{Float64}}() for t in algs)

        t_st_count = Dict{Float64, Int}()
        files = readdir(dirname)

        filter!(files) do fn
            isfile(joinpath(dirname, fn)) || return false
            startswith(fn, "stats_") && return false
            startswith(fn, "overlaps_") && return false
            startswith(fn, "Cs_") && return false
            ismatch(r"^output_(...)_sx\d+_s\d+\.txt$", fn) || (warn("skipping file $fn"); return false)
            tag = Symbol(match(r"output_(...)_sx\d+_s\d+\.txt", fn).captures[1])
            @assert tag ∈ algs # XXX
            return true
        end

        numfiles = length(files)
        fileind = 0
        filestr = ""
        for fn in files
            fileind += 1
            print("\r", " "^length(filestr), "\r")
            filestr = "analyzing file $fileind / $numfiles"
            print(filestr)

            tag = Symbol(match(r"output_(...)_sx\d+_s\d+\.txt", fn).captures[1])
            @assert tag ∈ algs # XXX

            # final_t = 0.0
            # final_gt = 1.0
            # if tag == :wtm
            #     ratio = open(joinpath(dirname, fn)) do f
            #         for l in eachline(f)
            #             eof(f) || continue
            #             sl = split(l)
            #             final_t = parse(Float64, sl[4])
            #             final_gt = parse(Float64, sl[1])
            #         end
            #         return final_t / final_gt
            #     end
            #     #info("ratio = $ratio")
            # else
            #     ratio = 1.0
            # end

            Ed_tag = Edict[tag]
            open(joinpath(dirname, fn)) do f
                max_file_t = 0.0
                Ed_tmp = Dict{Float64,Vector{Float64}}()
                for l in eachline(f)
                    ismatch(r"^\s*#", l) && continue
                    sl = split(l)
                    @assert length(sl) ≥ 4
                    E = parse(Float64, sl[3])
                    #t = tag == :wtm ? parse(Float64, sl[1]) * ratio : parse(Float64, sl[4])
                    t = parse(Float64, sl[4])

                    t_st = (t ÷ step) * step + step / 2
                    max_file_t = max(max_file_t, t_st)

                    Ed_tvec = get!(Ed_tmp, t_st, Float64[])
                    push!(Ed_tvec, E)
                end
                for t_st in sort!(collect(keys(Ed_tmp)))
                    Ed_tvec = Ed_tmp[t_st]
                    Ed_vec = get!(Ed_tag, t_st, Float64[])

                    Em = mean(Ed_tvec)
                    push!(Ed_vec, Em)

                    t_st_count[t_st] = get!(t_st_count, t_st, 0) + 1
                end
            end
        end
        println()

        #close(tf)

        max_t_st_count = maximum(values(t_st_count))
        all_t_st = sort!(collect(keys(t_st_count)))

        L = findfirst(t_st->t_st_count[t_st] < 1.0 * max_t_st_count, sort!(collect(keys(t_st_count))))
        L == 0 && (L = length(all_t_st))

        all_t_st = all_t_st[1:L]

        #all_t_st = filter!(t->t≤min_t, sort!(union([collect(keys(Ed_tag)) for Ed_tag in values(Edict)]...)))
        L = length(all_t_st)

        hasnans = false

        mdict = Dict{Symbol, Vector{Float64}}(tag => zeros(L) for tag in algs)
        zdict = Dict{Symbol, Vector{Float64}}(tag => zeros(L) for tag in algs)

        for tag in algs
            Ed_tag = Edict[tag]
            md_tag = mdict[tag]
            zd_tag = zdict[tag]

            for i = 1:L
                t_st = all_t_st[i]
                if !haskey(Ed_tag, t_st)
                    hasnans = true
                    md_tag[i] = NaN
                    zd_tag[i] = NaN
                    continue
                end

                md_tag[i] = mean(Ed_tag[t_st])
                zd_tag[i] = std(Ed_tag[t_st])
            end
        end

        hasnans && warn("empty bins found (increase step?)")

        for i = 1:L
            t_st = all_t_st[i]
            println(sf, t_st, " ", join((string(mdict[tag][i], " ", zdict[tag][i]) for tag in algs), " "))
        end
    end
end

# function stats_accrate_vs_β(dirname::AbstractString)
#     isdir(dirname) || error("directory $dirname not found")
#
#     tags = [:met, :rrr]
#
#     fname = joinpath(dirname, "stats_accrate_vs_beta.txt")
#
#     open(fname, "w") do sf
#         println(sf, "#β ", join((string(tag, " std") for tag in tags), " "))
#         sdict = Dict{Symbol, Dict{Float64, Vector{Float64}}}(t => Dict{Float64, Vector{Float64}}() for t in tags)
#
#         allβ = Set{Float64}()
#         files = readdir(dirname)
#         for fn in files
#             isfile(joinpath(dirname, fn)) || continue
#             startswith(fn, "stats_") && continue
#             ismatch(r"^output_(...)_s\d+\.txt$", fn) || (warn("skipping file $fn"); continue)
#             tag = Symbol(match(r"output_(...)_s\d+\.txt", fn).captures[1])
#             @assert tag ∈ tags
#             sd_tag = sdict[tag]
#             open(joinpath(dirname, fn)) do f
#                 β0 = -Inf
#                 it0 = 0
#                 it1 = 0
#                 acc1 = 0
#                 for l in eachline(f)
#                     ismatch(r"^\s*#", l) && continue
#                     sl = split(l)
#                     @assert length(sl) ≥ 5
#                     it = parse(Int, sl[1])
#                     acc = parse(Int, sl[2])
#                     E = parse(Float64, sl[3])
#                     t = parse(Float64, sl[4])
#                     β = parse(Float64, sl[5])
#
#                     if β ≠ β0
#                         push!(allβ, β)
#                         if β0 ≠ -Inf
#                             sd_tag_β = get!(sd_tag, β0) do
#                                 Float64[]
#                             end
#                             push!(sd_tag_β, acc1 / (it1 - it0))
#                         end
#                         it0 = it1
#                         β0 = β
#                     end
#                     it1 = it
#                     acc1 = acc
#                 end
#                 if β0 ≠ -Inf
#                     sd_tag_β = get!(sd_tag, β0) do
#                         Float64[]
#                     end
#                     push!(sd_tag_β, acc1 / (it1 - it0))
#                 end
#             end
#         end
#
#         vallβ = sort!(collect(allβ))
#         L = length(vallβ)
#
#         mdict = Dict{Symbol, Vector{Float64}}(tag => zeros(L) for tag in tags)
#         zdict = Dict{Symbol, Vector{Float64}}(tag => zeros(L) for tag in tags)
#
#         for tag in tags
#             sd_tag = sdict[tag]
#             md_tag = mdict[tag]
#             zd_tag = zdict[tag]
#
#             for i = 1:L
#                 β = vallβ[i]
#                 if !haskey(sd_tag, β)
#                     md_tag[i] = NaN
#                     zd_tag[i] = NaN
#                     continue
#                 end
#                 md_tag[i] = mean(sd_tag[β])
#                 zd_tag[i] = std(sd_tag[β])
#             end
#         end
#
#         for i = 1:L
#             β = vallβ[i]
#             println(sf, β, " ", join((string(mdict[tag][i], " ", zdict[tag][i]) for tag in tags), " "))
#         end
#     end
# end

# function stats_accrate_vs_time(dirname::AbstractString; step::Float64 = 10.0)
#     isdir(dirname) || error("directory $dirname not found")
#
#     tags = [:met, :rrr]
#
#     fname = joinpath(dirname, "stats_accrate_vs_time.txt")
#
#     open(fname, "w") do sf
#         println(sf, "#time ", join((string(tag, " std") for tag in tags), " "))
#         sdict = Dict{Symbol, Dict{Float64, Vector{Float64}}}(t => Dict{Float64, Vector{Float64}}() for t in tags)
#
#         min_t = Inf
#         files = readdir(dirname)
#         for fn in files
#             isfile(joinpath(dirname, fn)) || continue
#             startswith(fn, "stats_") && continue
#             ismatch(r"^output_(...)_s\d+\.txt$", fn) || (warn("skipping file $fn"); continue)
#             tag = Symbol(match(r"output_(...)_s\d+\.txt", fn).captures[1])
#             @assert tag ∈ tags
#             sd_tag = sdict[tag]
#             open(joinpath(dirname, fn)) do f
#                 max_file_t = 0.0
#                 it0 = 0
#                 acc0 = 0
#                 for l in eachline(f)
#                     ismatch(r"^\s*#", l) && continue
#                     sl = split(l)
#                     @assert length(sl) ≥ 4
#                     it = parse(Int, sl[1])
#                     acc = parse(Int, sl[2])
#                     E = parse(Float64, sl[3])
#                     t = parse(Float64, sl[4])
#
#                     t_st = (t ÷ step) * step + step / 2
#                     max_file_t = max(max_file_t, t_st)
#
#                     acc_r = (acc - acc0) / (it - it0)
#
#                     acc0, it0 = acc, it
#
#                     sd_vec = get!(sd_tag, t_st) do
#                         Float64[]
#                     end
#                     push!(sd_vec, acc_r)
#                 end
#                 min_t = min(min_t, max_file_t)
#             end
#         end
#
#         all_t_st = filter!(t->t≤min_t, sort!(union([collect(keys(sd_tag)) for sd_tag in values(sdict)]...)))
#         L = length(all_t_st)
#
#         hasnans = false
#
#         mdict = Dict{Symbol, Vector{Float64}}(tag => zeros(L) for tag in tags)
#         zdict = Dict{Symbol, Vector{Float64}}(tag => zeros(L) for tag in tags)
#
#         for tag in tags
#             sd_tag = sdict[tag]
#             md_tag = mdict[tag]
#             zd_tag = zdict[tag]
#
#             for i = 1:L
#                 t_st = all_t_st[i]
#                 if !haskey(sd_tag, t_st)
#                     hasnans = true
#                     md_tag[i] = NaN
#                     zd_tag[i] = NaN
#                     continue
#                 end
#                 md_tag[i] = mean(sd_tag[t_st])
#                 zd_tag[i] = std(sd_tag[t_st])
#             end
#         end
#
#         hasnans && warn("empty bins found (increase step?)")
#
#         for i = 1:L
#             t_st = all_t_st[i]
#             println(sf, t_st, " ", join((string(mdict[tag][i], " ", zdict[tag][i]) for tag in tags), " "))
#         end
#     end
# end

end # module
