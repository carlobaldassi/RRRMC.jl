module RRRMCTest

using RRRMC
using Base.Test

function gen_timeout_hook(t = 1.0)
    t += time()
    return (args...) -> (time() ≤ t)
end

function checkenergy_hook(it, X, C, acc, E)
    @test ≈(E, RRRMC.energy(X, C), atol=1e-12)
    # @test E ≈ RRRMC.energy(X, C) atol=1e-12 # change to this when julia v0.5 support is dropped
    return true
end

function checkenergy_hook_EO(it, X, C, E, Emin)
    @test ≈(E, RRRMC.energy(X, C), atol=1e-12)
    # @test E ≈ RRRMC.energy(X, C) atol=1e-12 # change to this when julia v0.5 support is dropped
    return true
end

function test()

    srand(8426732438942) # get reproducible results...

    graphs = [
        RRRMC.GraphTwoSpin(),
        RRRMC.GraphThreeSpin(),
        RRRMC.GraphFields(10),
        RRRMC.GraphFields(10, (-1,1)),
        RRRMC.GraphFields(10, (-1.5,0.5)),
        RRRMC.GraphFieldsNormalDiscretized(10, (-1,0,1)),
        RRRMC.GraphFieldsNormalDiscretized(10, (-1.5,0.0,1.5)),
        RRRMC.GraphIsing1D(10),
        RRRMC.GraphPSpin3(39, 5),
        RRRMC.GraphRRG(10, 3),
        RRRMC.GraphRRG(10, 3, (-1,0,1)),
        RRRMC.GraphRRG(10, 3, (-1.0,0.0,1.0)),
        RRRMC.GraphRRG(10, 3, (-1.0,0,1)),
        RRRMC.GraphRRG(10, 3, (-1//1,0//1,1//1)),
        RRRMC.GraphRRGNormalDiscretized(10, 3, (-1,0,1)),
        RRRMC.GraphRRGNormalDiscretized(10, 3, (-1.0,0.0,1.0)),
        RRRMC.GraphRRGNormalDiscretized(10, 3, (-1.0,0,1)),
        RRRMC.GraphRRGNormalDiscretized(10, 3, (-1//1,0//1,1//1)),
        RRRMC.GraphRRGNormal(10, 3),
        RRRMC.GraphEA(2, 3),
        RRRMC.GraphEA(2, 3, (-1,0,1)),
        RRRMC.GraphEA(2, 3, (-1.0,0.0,1.0)),
        RRRMC.GraphEA(2, 3, (-1.0,0,1)),
        RRRMC.GraphEA(2, 3, (-1//1,0//1,1//1)),
        RRRMC.GraphEANormalDiscretized(2, 3, (-1,0,1)),
        RRRMC.GraphEANormalDiscretized(2, 3, (-1.0,0.0,1.0)),
        RRRMC.GraphEANormalDiscretized(2, 3, (-1.0,0,1)),
        RRRMC.GraphEANormalDiscretized(2, 3, (-1//1,0//1,1//1)),
        RRRMC.GraphEANormal(2, 3),
        RRRMC.GraphEA(3, 2),
        RRRMC.GraphEA(3, 2, (-1,0,1)),
        RRRMC.GraphEA(3, 2, (-1.0,0.0,1.0)),
        RRRMC.GraphEA(3, 2, (-1.0,0,1)),
        RRRMC.GraphEA(3, 2, (-1//1,0//1,1//1)),
        RRRMC.GraphEANormalDiscretized(3, 2, (-1,0,1)),
        RRRMC.GraphEANormalDiscretized(3, 2, (-1.0,0.0,1.0)),
        RRRMC.GraphEANormalDiscretized(3, 2, (-1.0,0,1)),
        RRRMC.GraphEANormalDiscretized(3, 2, (-1//1,0//1,1//1)),
        RRRMC.GraphEANormal(3, 2),
        RRRMC.GraphSK(10),
        RRRMC.GraphSKNormal(10),
        RRRMC.GraphSAT(10, 3, 4.2),
        RRRMC.GraphPercLinear(101, 30),
        RRRMC.GraphPercStep(101, 30),
        RRRMC.GraphCommStep(21, 5, 30),
        RRRMC.GraphQuant(10, 8, 0.5, 2.0, RRRMC.GraphEmpty, 10),
        RRRMC.GraphQuant(10, 8, 0.5, 2.0, RRRMC.GraphSK, RRRMC.SK.gen_J(10)),
        RRRMC.GraphQuant(10, 8, 0.5, 2.0, RRRMC.GraphSKNormal, RRRMC.SK.gen_J_gauss(10)),
        RRRMC.GraphQuant(3, 8, 0.5, 2.0, RRRMC.GraphThreeSpin),
        RRRMC.GraphQPercLinearT(101, 30, 5, 0.5, 2.0),
        RRRMC.GraphQPercStepT(101, 30, 5, 0.5, 2.0),
        RRRMC.GraphQCommStepT(21, 5, 30, 5, 0.5, 2.0),
        RRRMC.GraphRobustEnsemble(10, 8, 1.5, 2.0, RRRMC.GraphEmpty, 10),
        RRRMC.GraphRobustEnsemble(10, 8, 1.5, 2.0, RRRMC.GraphSK, RRRMC.SK.gen_J(10)),
        RRRMC.GraphRobustEnsemble(10, 8, 1.5, 2.0, RRRMC.GraphSKNormal, RRRMC.SK.gen_J_gauss(10)),
        RRRMC.GraphRobustEnsemble(3, 8, 1.5, 2.0, RRRMC.GraphThreeSpin),
        RRRMC.GraphPercLinearRE(101, 30, 5, 0.5, 2.0),
        RRRMC.GraphPercStepRE(101, 30, 5, 0.5, 2.0),
        RRRMC.GraphCommStepRE(21, 5, 30, 5, 0.5, 2.0),
        RRRMC.GraphLocalEntropy(10, 8, 1.5, 2.0, RRRMC.GraphEmpty, 10),
        RRRMC.GraphLocalEntropy(10, 8, 1.5, 2.0, RRRMC.GraphSKNormal, RRRMC.SK.gen_J_gauss(10)),
        RRRMC.GraphLocalEntropy(3, 8, 1.5, 2.0, RRRMC.GraphThreeSpin),
        RRRMC.GraphRobustEnsemble(20, 4, 1.5, 2.0, RRRMC.GraphQuant, 5, 4, 0.5, 2.0, RRRMC.GraphSK, RRRMC.SK.gen_J(5)),
        RRRMC.GraphPercLinearLE(101, 30, 5, 0.5, 2.0),
        RRRMC.GraphPercStepLE(101, 30, 5, 0.5, 2.0),
        RRRMC.GraphCommStepLE(21, 5, 30, 5, 0.5, 2.0),
       ]

    β = 2.0
    iters = 10_000
    st = 100
    τ = 1.3
    samples = iters ÷ st
    quiet = true

    for X in graphs
        # @show X
        E, C = standardMC(X, β, iters, step=st, quiet=quiet)
        E, C = standardMC(X, β, iters, step=st, quiet=quiet, C0=C, hook=checkenergy_hook)
        E, C = standardMC(X, β, iters, step=st, quiet=quiet, C0=C, hook=gen_timeout_hook())

        E, C = bklMC(X, β, iters, step=st, quiet=quiet)
        E, C = bklMC(X, β, iters, step=st, quiet=quiet, C0=C, hook=checkenergy_hook)
        E, C = bklMC(X, β, iters, step=st, quiet=quiet, C0=C, hook=gen_timeout_hook())

        E, C = wtmMC(X, β, samples, step=Float64(st), quiet=quiet)
        E, C = wtmMC(X, β, samples, step=Float64(st), quiet=quiet, C0=C, hook=checkenergy_hook)
        E, C = wtmMC(X, β, samples, step=Float64(st), quiet=quiet, C0=C, hook=gen_timeout_hook())

        E, C = rrrMC(X, β, iters, step=st, quiet=quiet)
        E, C = rrrMC(X, β, iters, step=st, quiet=quiet, C0=C, hook=checkenergy_hook)
        E, C = rrrMC(X, β, iters, step=st, quiet=quiet, C0=C, hook=gen_timeout_hook())
        E, C = rrrMC(X, β, iters, step=st, quiet=quiet, staged_thr=0.0, hook=checkenergy_hook)
        E, C = rrrMC(X, β, iters, step=st, quiet=quiet, C0=C, staged_thr=0.0)
        E, C = rrrMC(X, β, iters, step=st, quiet=quiet, staged_thr=1.0, hook=checkenergy_hook)
        E, C = rrrMC(X, β, iters, step=st, quiet=quiet, C0=C, staged_thr=1.0)

        C, Emin, Cmin, itmin = extremal_opt(X, τ, iters, step=st, quiet=quiet)
        C, Emin, Cmin, itmin = extremal_opt(X, τ, iters, step=st, quiet=quiet, C0=C)
        C, Emin, Cmin, itmin = extremal_opt(X, τ, iters, step=st, quiet=quiet, hook=checkenergy_hook_EO)
        C, Emin, Cmin, itmin = extremal_opt(X, τ, iters, step=st, quiet=quiet, hook=gen_timeout_hook())

        if isa(X, RRRMC.DoubleGraph)
            X0 = RRRMC.inner_graph(X)
            E, C = bklMC(X0, β, iters, step=st, quiet=quiet)
            E, C = bklMC(X0, β, iters, step=st, quiet=quiet, C0=C, hook=checkenergy_hook)

            E, C = rrrMC(X, β, iters, step=st, quiet=quiet)
            E, C = rrrMC(X, β, iters, step=st, quiet=quiet, C0=C, hook=checkenergy_hook)
            E, C = rrrMC(X, β, iters, step=st, quiet=quiet, C0=C, hook=gen_timeout_hook())
            E, C = rrrMC(X, β, iters, step=st, quiet=quiet, staged_thr=0.0, hook=checkenergy_hook)
            E, C = rrrMC(X, β, iters, step=st, quiet=quiet, C0=C, staged_thr=0.0)
            E, C = rrrMC(X, β, iters, step=st, quiet=quiet, staged_thr=1.0, hook=checkenergy_hook)
            E, C = rrrMC(X, β, iters, step=st, quiet=quiet, C0=C, staged_thr=1.0)


            E, C = rrrMC(X0, β, iters, step=st, quiet=quiet)
            E, C = rrrMC(X0, β, iters, step=st, quiet=quiet, C0=C, hook=checkenergy_hook)
            E, C = rrrMC(X0, β, iters, step=st, quiet=quiet, C0=C, hook=gen_timeout_hook())
            E, C = rrrMC(X0, β, iters, step=st, quiet=quiet, staged_thr=0.0, hook=checkenergy_hook)
            E, C = rrrMC(X0, β, iters, step=st, quiet=quiet, C0=C, staged_thr=0.0)
            E, C = rrrMC(X0, β, iters, step=st, quiet=quiet, staged_thr=1.0, hook=checkenergy_hook)
            E, C = rrrMC(X0, β, iters, step=st, quiet=quiet, C0=C, staged_thr=1.0)

            C, Emin, Cmin, itmin = extremal_opt(X0, τ, iters, step=st, quiet=quiet)
            C, Emin, Cmin, itmin = extremal_opt(X0, τ, iters, step=st, quiet=quiet, C0=C)
            C, Emin, Cmin, itmin = extremal_opt(X0, τ, iters, step=st, quiet=quiet, hook=checkenergy_hook_EO)
            C, Emin, Cmin, itmin = extremal_opt(X0, τ, iters, step=st, quiet=quiet, hook=gen_timeout_hook())
        end
    end
end

test()

end # module
