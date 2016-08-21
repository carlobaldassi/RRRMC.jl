module RRRMCTest

using RRRMC
using Base.Test

function test()

    graphs = [
        RRRMC.GraphTwoSpin(),
        RRRMC.GraphThreeSpin(),
        RRRMC.GraphFields(10),
        RRRMC.GraphFields(10, (-1,1)),
        RRRMC.GraphFields(10, (-1.5,0.5)),
        RRRMC.GraphFieldsCont(10, (-1,0,1)),
        RRRMC.GraphFieldsCont(10, (-1.5,0.0,1.5)),
        RRRMC.GraphIsing1D(10),
        RRRMC.GraphPSpin3(39, 5),
        RRRMC.GraphRRG(10, 3),
        RRRMC.GraphRRG(10, 3, (-1,0,1)),
        RRRMC.GraphRRG(10, 3, (-1.0,0.0,1.0)),
        RRRMC.GraphRRG(10, 3, (-1.0,0,1)),
        RRRMC.GraphRRG(10, 3, (-1//1,0//1,1//1)),
        RRRMC.GraphRRGCont(10, 3, (-1,0,1)),
        RRRMC.GraphRRGCont(10, 3, (-1.0,0.0,1.0)),
        RRRMC.GraphRRGCont(10, 3, (-1.0,0,1)),
        RRRMC.GraphRRGCont(10, 3, (-1//1,0//1,1//1)),
        RRRMC.GraphRRGContSimple(10, 3),
        RRRMC.GraphEA(2, 3),
        RRRMC.GraphEA(2, 3, (-1,0,1)),
        RRRMC.GraphEA(2, 3, (-1.0,0.0,1.0)),
        RRRMC.GraphEA(2, 3, (-1.0,0,1)),
        RRRMC.GraphEA(2, 3, (-1//1,0//1,1//1)),
        RRRMC.GraphEACont(2, 3, (-1,0,1)),
        RRRMC.GraphEACont(2, 3, (-1.0,0.0,1.0)),
        RRRMC.GraphEACont(2, 3, (-1.0,0,1)),
        RRRMC.GraphEACont(2, 3, (-1//1,0//1,1//1)),
        RRRMC.GraphEAContSimple(2, 3),
        RRRMC.GraphEA(3, 2),
        RRRMC.GraphEA(3, 2, (-1,0,1)),
        RRRMC.GraphEA(3, 2, (-1.0,0.0,1.0)),
        RRRMC.GraphEA(3, 2, (-1.0,0,1)),
        RRRMC.GraphEA(3, 2, (-1//1,0//1,1//1)),
        RRRMC.GraphEACont(3, 2, (-1,0,1)),
        RRRMC.GraphEACont(3, 2, (-1.0,0.0,1.0)),
        RRRMC.GraphEACont(3, 2, (-1.0,0,1)),
        RRRMC.GraphEACont(3, 2, (-1//1,0//1,1//1)),
        RRRMC.GraphEAContSimple(3, 2),
        RRRMC.GraphQ0T(10, 8, 0.5, 2.0),
        RRRMC.GraphQIsingT(10, 8, 0.5, 2.0)
       ]

    β = 2.0
    iters = 10_000
    st = 100

    for X in graphs
        E, C = standardMC(X, β, iters, step=st)
        E, C = standardMC(X, β, iters, step=st, C0=C)

        if isa(X, RRRMC.DiscrGraph)
            E, C = bklMC(X, β, iters, step=st)
            E, C = bklMC(X, β, iters, step=st, C0=C)

            E, C = rrrMC(X, β, iters, step=st)
            E, C = rrrMC(X, β, iters, step=st, C0=C)
            E, C = rrrMC(X, β, iters, step=st, staged_thr=0.0)
            E, C = rrrMC(X, β, iters, step=st, C0=C, staged_thr=0.0)
            E, C = rrrMC(X, β, iters, step=st, staged_thr=1.0)
            E, C = rrrMC(X, β, iters, step=st, C0=C, staged_thr=1.0)
        elseif isa(X, RRRMC.DoubleGraph)
            X0 = RRRMC.discr_graph(X)
            E, C = bklMC(X0, β, iters, step=st)
            E, C = bklMC(X0, β, iters, step=st, C0=C)

            E, C = rrrMC(X, β, iters, step=st)
            E, C = rrrMC(X, β, iters, step=st, C0=C)
            E, C = rrrMC(X, β, iters, step=st, staged_thr=0.0)
            E, C = rrrMC(X, β, iters, step=st, C0=C, staged_thr=0.0)
            E, C = rrrMC(X, β, iters, step=st, staged_thr=1.0)
            E, C = rrrMC(X, β, iters, step=st, C0=C, staged_thr=1.0)


            E, C = rrrMC(X0, β, iters, step=st)
            E, C = rrrMC(X0, β, iters, step=st, C0=C)
            E, C = rrrMC(X0, β, iters, step=st, staged_thr=0.0)
            E, C = rrrMC(X0, β, iters, step=st, C0=C, staged_thr=0.0)
            E, C = rrrMC(X0, β, iters, step=st, staged_thr=1.0)
            E, C = rrrMC(X0, β, iters, step=st, C0=C, staged_thr=1.0)
        end
    end
end

test()

end # module
