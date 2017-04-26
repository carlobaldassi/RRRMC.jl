var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#RRRMC.jl-documentation-1",
    "page": "Home",
    "title": "RRRMC.jl documentation",
    "category": "section",
    "text": "CurrentModule = RRRMCThis code implements the Reduced-Rejection-Rate (RRR) Monte Carlo method for Ising spin models described in the paper \"A method to reduce the rejection rate in Monte Carlo Markov Chains\" by C. Baldassi, J. Stat. Mech. Theor. Exp., (2017) 3, 033301 [doi:10.1088/1742-5468/aa5335]https://doi.org/10.1088/1742-5468/aa5335.It also provides:a standard Metropolis-Hastings sampler\na generalized implementation of the BKL method described in the paper \"A new algorithm for Monte Carlo simulation of Ising spin systems\" by A.B. Bortz, M.H. Kalos and J.L. Lebowitz. The generalization consists in not requiring that the energy shifts are discrete.\nan implementation of the Waiting time method described in the paper \"Faster Monte Carlo simulations at low temperatures. The waiting time method\" by J. Dall and P. Sibani.The code is written in Julia, and tested against Julia 0.5 and current 0.6-dev on Linux, OS X, and Windows."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "To install the module, use this command from within Julia:julia> Pkg.clone(\"https://github.com/carlobaldassi/RRRMC.jl\")Dependencies will be installed automatically."
},

{
    "location": "index.html#Usage-1",
    "page": "Home",
    "title": "Usage",
    "category": "section",
    "text": "The module is loaded as any other Julia module:julia> using RRRMCThe module provides four functions which implement Monte Carlo Markov Chain algorithms on Ising spin models:standardMC: a standard Metropolis-Hastings sampler\nrrrMC: the reduced-rejection-rate (RRR) method\nbklMC: the Bortz-Kalos-Lebowitz (BKL) method\nwtmMC: the waiting-time method (WTM)The interface for these four algorithms is documented in the Sampling algorithms page, and it is essentially identical: they take as arguments a graph, an inverse temperature parameter β, and the number of Monte Carlo iterations to perform (or, for wtmMC, of samples to collect). However, the sampling methodology changes based on the type of model, see the Graph types page.These functions allow accessing the internal state during the iteration at regular intervals, via the hook keyword argument. They also return the final configuration of the system, which is stored in an object of type Config.The code comes with some built-in graphs, but provides an interface to write user-defined models.note: Note\nThe four sampling functions are the only names exported by the module; all other function and types must be qualified with the RRRMC module name."
},

{
    "location": "index.html#Manual-1",
    "page": "Home",
    "title": "Manual",
    "category": "section",
    "text": "Pages = [\n    \"algorithms.md\",\n    \"graph-types.md\",\n    \"graphs-builtin.md\",\n    \"interface.md\",\n    ]\nDepth = 3"
},

{
    "location": "algorithms.html#",
    "page": "Sampling algorithms",
    "title": "Sampling algorithms",
    "category": "page",
    "text": ""
},

{
    "location": "algorithms.html#RRRMC.standardMC",
    "page": "Sampling algorithms",
    "title": "RRRMC.standardMC",
    "category": "Function",
    "text": "standardMC(X::AbstractGraph, β::Real, iters::Integer; keywords...)\n\nRuns iters iterations of a standard Metropolis Monte Carlo algorithm on the given Ising spin model X, at inverse temperature β. Each spin flip attempt counts as an iteration.\n\nReturns two objects: a vector of energies and the last configuration (see Config).\n\nPossible keyord arguments are:\n\nstep: the interval of iterations with which to collect energies (for the returned result) and to call the hook function (see below). The default is 1, which is good for debugging but otherwise generally a bad idea.\nseed: the random seed. The default is some arbitrary number.\nC0: the initial configuration. The default is nothing, in which case it is initialized at random. Otherwise it can be a Config object. Passing the result of a previous run can be useful e.g. when implementing a simulated annealing protocol, or if the system has not equilibrated yet.\nhook: a function to be executed after every step number of iterations (see above). It must take five arguments: the current iteration, the graph X, the current configuration, the number of accepted moves so far, and the current energy. Useful to collect data other than the energy, write to files ecc; you'd probably want to use a closure, see the example below. The return value must be a Bool: return false to interrupt the simulation, true otherwise. The default does nothing and just returns true.\n\nBasic example:\n\njulia> srand(76543); X = RRRMC.GraphPSpin3(3999, 5); β = 1.0;\njulia> Es, C = standardMC(X, β, 100_000, step = 1_000);\n\nExample of using the hook for collecting samples as the columns of a BitMatrix:\n\njulia> iters = 100_000; step = 1_000; l = iters ÷ step; N = RRRMC.getN(X);\njulia> Cs = BitArray(N, l); hook = (it, X, C, acc, E) -> (Cs[:,it÷step]=C.s; true);\njulia> Es, C = standardMC(X, β, iters, step = step, hook = hook);\n\n\n\n"
},

{
    "location": "algorithms.html#RRRMC.rrrMC",
    "page": "Sampling algorithms",
    "title": "RRRMC.rrrMC",
    "category": "Function",
    "text": "rrrMC(X::AbstractGraph, β::Real, iters::Integer; keywords...)\n\nSame as standardMC, but uses the reduced-rejection-rate method. Each iteration takes more time, but has a higher chance of being accepted, so fewer iterations overall should be needed normally. Whether this trade-off is convenient depends on the parameters and the details of the model. This function has specialized versions for DiscrGraph and DoubleGraph models.\n\nThe return values and the keyword arguments are the same as standardMC, see the usage examples for that function.\n\n\n\n"
},

{
    "location": "algorithms.html#RRRMC.bklMC",
    "page": "Sampling algorithms",
    "title": "RRRMC.bklMC",
    "category": "Function",
    "text": "bklMC(X::AbstractGraph, β::Real, iters::Integer; keywords...)\n\nSame as standardMC, but uses the rejection-free method by Bortz, Kalos and Lebowitz. Each step takes more time, but rejected moves are almost free, since they are skipped entirely. This function has a specialized version for DiscrGraph models.\n\nThe return values and the keyword arguments are the same as standardMC, see the usage examples for that function.\n\nNote that the number of iterations includes the rejected moves. This makes the results directly comparable with those of standardMC. It also means that increasing β at fixed iters will result in fewer steps being actually computed.\n\n\n\n"
},

{
    "location": "algorithms.html#RRRMC.wtmMC",
    "page": "Sampling algorithms",
    "title": "RRRMC.wtmMC",
    "category": "Function",
    "text": "wtmMC(X::AbstractGraph, β::Real, samples::Integer; keywords...)\n\nSame as standardMC, but uses the rejection-free waiting-time method by Dall and Sibani. It is similar to bklMC.\n\nThe return values and the keyword arguments are almost the same as standardMC, see the usage examples for that function. However, the waiting time method uses an internal \"global time\" variable, which takes the place of the \"iterations\" counter of standardMC and of bklMC. Thus, this function has two differences with respect to the other samplers in the module:\n\nthe function takes a samples integer argument, with the maximum number of samples which will be collected.\nthe step keyword argument is of type Float64 instead of integer (default value = 1.0; you'll probably want to change this). The step is measured in terms of the global time variable and is scaled with the size of the problem N.\n\nThe total number of samples actually collected can still be less than samples if the hook function from the keyword arguments returns false earlier.\n\n\n\n"
},

{
    "location": "algorithms.html#RRRMC.Interface.Config",
    "page": "Sampling algorithms",
    "title": "RRRMC.Interface.Config",
    "category": "Type",
    "text": "Config(N::Integer)\n\nThe object storing the configuration for an Ising model. Although the spin values are _i  -11, internally they are stored in a BitArray, in the type field s, so that to obtain the real value one needs to perform the transformation _i = 2s_i - 1.\n\n\n\n"
},

{
    "location": "algorithms.html#algorithms-1",
    "page": "Sampling algorithms",
    "title": "Sampling algorithms",
    "category": "section",
    "text": "CurrentModule = RRRMCstandardMCrrrMCbklMCwtmMCConfig"
},

{
    "location": "graph-types.html#",
    "page": "Graph types",
    "title": "Graph types",
    "category": "page",
    "text": ""
},

{
    "location": "graph-types.html#RRRMC.Interface.AbstractGraph",
    "page": "Graph types",
    "title": "RRRMC.Interface.AbstractGraph",
    "category": "Type",
    "text": "AbstractGraph{ET<:Real}\n\nAn abstract type representing an Ising spin model. The ET parameter is the type returned by the energy and delta_energy functions.\n\nSee also SimpleGraph, DiscrGraph, SingleGraph and DoubleGraph.\n\n\n\n"
},

{
    "location": "graph-types.html#RRRMC.Interface.SimpleGraph",
    "page": "Graph types",
    "title": "RRRMC.Interface.SimpleGraph",
    "category": "Type",
    "text": "SimpleGraph{ET} <: AbstractGraph{ET}\n\nAn abstract type representing a generic graph.\n\nThe ET parameter is the type returned by energy and delta_energy.\n\n\n\n"
},

{
    "location": "graph-types.html#RRRMC.Interface.DiscrGraph",
    "page": "Graph types",
    "title": "RRRMC.Interface.DiscrGraph",
    "category": "Type",
    "text": "DiscrGraph{ET} <: AbstractGraph{ET}\n\nAn abstract type representing a graph in which the delta_energy values produced when flipping a spin belong to a finite discrete set, and thus can be sampled more efficiently with rrrMC or bklMC if the set is small.\n\nThe ET parameter is the type returned by energy and delta_energy.\n\nSee also allΔE.\n\n\n\n"
},

{
    "location": "graph-types.html#RRRMC.Interface.SingleGraph",
    "page": "Graph types",
    "title": "RRRMC.Interface.SingleGraph",
    "category": "Constant",
    "text": "SingleGraph{ET}\n\nA type alias representing either a SimpleGraph or a DiscrGraph{ET}. See also DoubleGraph.\n\n\n\n"
},

{
    "location": "graph-types.html#RRRMC.Interface.DoubleGraph",
    "page": "Graph types",
    "title": "RRRMC.Interface.DoubleGraph",
    "category": "Type",
    "text": "DoubleGraph{GT<:SingleGraph,ET} <: AbstractGraph{ET}\n\nAn abstract type representing a graph in which the energy is the sum of two contributions, one of which is encoded in a graph of type GT (see SingleGraph). This allows rrrMC to sample values more efficiently.\n\nThe ET parameter is the type returned by the energy and delta_energy functions. Note that it can be different from the energy type of the internal GT object (e.g., one can have a DoubleGraph{DiscrGraph{Int},Float64} object).\n\nNote: When you declare a type as subtype of this, GT should not be the concrete type of the inner graph, but either SimpleGraph{T} or DiscrGraph{T} for some T.\n\nSee also inner_graph, delta_energy_residual and update_cache_residual!.\n\n\n\n"
},

{
    "location": "graph-types.html#graphtype-1",
    "page": "Graph types",
    "title": "Graph types",
    "category": "section",
    "text": "CurrentModule = RRRMCAll graphs which can be used with the sampling algorithms belong to a type hierarchy. At the top of the hierarchy, there is AbstractGraph:AbstractGraphThere are currently three abstract subclasses and an alias, which determine how the sampling algorithms can be used:SimpleGraphDiscrGraphSingleGraphDoubleGraph"
},

{
    "location": "graphs-builtin.html#",
    "page": "Built-in graphs",
    "title": "Built-in graphs",
    "category": "page",
    "text": ""
},

{
    "location": "graphs-builtin.html#builtin-1",
    "page": "Built-in graphs",
    "title": "Built-in graphs",
    "category": "section",
    "text": "CurrentModule = RRRMCFollowing is the list of the graph models which are provided with the module. After loading the RRRMC module, they can be constructed like in this example:julia> X = RRRMC.GraphRRG(10, 3)Note that for models which involve randomness in the constructor you may want to set the random seed with srand before calling the constructor, for reproducibility purposes."
},

{
    "location": "graphs-builtin.html#Basic-spin-glass-models-1",
    "page": "Built-in graphs",
    "title": "Basic spin glass models",
    "category": "section",
    "text": ""
},

{
    "location": "graphs-builtin.html#RRRMC.RRG.GraphRRG",
    "page": "Built-in graphs",
    "title": "RRRMC.RRG.GraphRRG",
    "category": "Type",
    "text": "GraphRRG(N::Integer, K::Integer, LEV = (-1,1)) <: DiscrGraph\n\nA DiscGraph implementing a random regular graph with N spins and connectivity K. Note: N*K must be even. Also, the graph generator uses the pairing model method by Bollobás, with a cutoff on the number of restarts, and thus it may occasionally fail if K is large.\n\nThe interactions are extracted at random from LEV, which must be a Tuple of Reals. No external fields.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#RRRMC.RRG.GraphRRGNormal",
    "page": "Built-in graphs",
    "title": "RRRMC.RRG.GraphRRGNormal",
    "category": "Type",
    "text": "GraphRRGNormal(N::Integer, K::Integer) <: SimpleGraph{Float64}\n\nA SimpleGraph implementing a random regular graph with N spins and connectivity K. Note: N*K must be even. Also, the graph generator uses the pairing model method by Bollobás, with a cutoff on the number of restarts, and thus it may occasionally fail if K is large.\n\nSame as GraphRRG, but the interactions are extracted from a normal distribution with unit variance.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#RRRMC.RRG.GraphRRGNormalDiscretized",
    "page": "Built-in graphs",
    "title": "RRRMC.RRG.GraphRRGNormalDiscretized",
    "category": "Type",
    "text": "GraphRRGNormalDiscretized(N::Integer, K::Integer, LEV) <: DoubleGraph{DiscrGraph,Float64}\n\nA DoubleGraph implementing a random regular graph with N spins and connectivity K. Note: N*K must be even. Also, the graph generator uses the pairing model method by Bollobás, with a cutoff on the number of restarts, and thus it may occasionally fail if K is large.\n\nThe interactions are extracted from a normal distribution with unit variance, and are then discretized using the values in LEV, which must be a Tuple of Reals. No external fields.\n\nSame as GraphRRGNormal, but it works differently when used with rrrMC.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#Random-regular-graphs-1",
    "page": "Built-in graphs",
    "title": "Random regular graphs",
    "category": "section",
    "text": "GraphRRGGraphRRGNormalGraphRRGNormalDiscretized"
},

{
    "location": "graphs-builtin.html#RRRMC.EA.GraphEA",
    "page": "Built-in graphs",
    "title": "RRRMC.EA.GraphEA",
    "category": "Type",
    "text": "GraphEA(L::Integer, D::Integer, LEV = (-1,1)) <: DiscrGraph\n\nAn Edwards-Anderson DiscrGraph: spins are arranged on a square lattice of size L in D dimensions (i.e. there are L^D total spins), with periodic boundary conditions.\n\nThe interactions are extracted at random from LEV, which must be a Tuple of Reals. No external fields.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#RRRMC.EA.GraphEANormal",
    "page": "Built-in graphs",
    "title": "RRRMC.EA.GraphEANormal",
    "category": "Type",
    "text": "GraphEANormal(L::Integer, D::Integer) <: SimpleGraph{Float64}\n\nAn Edwards-Anderson SimpleGraph: spins are arranged on a square lattice of size L in D dimensions (i.e. there are L^D total spins), with periodic boundary conditions.\n\nSame as GraphEA, but the interactions are extracted from a normal distribution with unit variance.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#RRRMC.EA.GraphEANormalDiscretized",
    "page": "Built-in graphs",
    "title": "RRRMC.EA.GraphEANormalDiscretized",
    "category": "Type",
    "text": "GraphEANormalDiscretized(L::Integer, D::Integer, LEV) <: DoubleGraph{DiscrGraph,Float64}\n\nAn Edwards-Anderson DoubleGraph: spins are arranged on a square lattice of size L in D dimensions (i.e. there are L^D total spins), with periodic boundary conditions.\n\nThe interactions are extracted at random from a normal distribution with unit variance, and are then discretized using the values in LEV, which must be a Tuple of Reals. No external fields.\n\nSame as GraphEANormal, but works differently when used with rrrMC.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#Edwards-Anderson-graphs-1",
    "page": "Built-in graphs",
    "title": "Edwards-Anderson graphs",
    "category": "section",
    "text": "GraphEAGraphEANormalGraphEANormalDiscretized"
},

{
    "location": "graphs-builtin.html#RRRMC.PSpin3.GraphPSpin3",
    "page": "Built-in graphs",
    "title": "RRRMC.PSpin3.GraphPSpin3",
    "category": "Type",
    "text": "GraphPSpin3(N::Integer, K::Integer) <: DiscrGraph{Int}\n\nA DiscrGraph implementing a p-spin regular graph with p=3. N is the number of spins, and must be divisible by 3; K is the connectivity. All interactions are set to J=1.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#p-spin-1",
    "page": "Built-in graphs",
    "title": "p-spin",
    "category": "section",
    "text": "GraphPSpin3"
},

{
    "location": "graphs-builtin.html#RRRMC.SK.GraphSK",
    "page": "Built-in graphs",
    "title": "RRRMC.SK.GraphSK",
    "category": "Type",
    "text": "GraphSK(N::Integer) <: SimpleGraph{Float64}\n\nA SimpleGraph implementing a Sherrington-Kirkpatrick fully-connected Ising model with N spins and random binary interactions (J  -1N1N) and no external fields.\n\nSame as GraphSKNormal, but with binary interactions.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#RRRMC.SK.GraphSKNormal",
    "page": "Built-in graphs",
    "title": "RRRMC.SK.GraphSKNormal",
    "category": "Type",
    "text": "GraphSK(N::Integer) <: SimpleGraph{Float64}\n\nA SimpleGraph implementing a Sherrington-Kirkpatrick fully-connected Ising model with N spins and random interactions extracted from a normal distribution with zero mean and 1N variance, and no external fields.\n\nSame as GraphSK, but with Gaussian interactions.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#Sherrington-Kirkpatrick-graphs-1",
    "page": "Built-in graphs",
    "title": "Sherrington-Kirkpatrick graphs",
    "category": "section",
    "text": "GraphSKGraphSKNormal"
},

{
    "location": "graphs-builtin.html#RRRMC.QT.GraphQuant",
    "page": "Built-in graphs",
    "title": "RRRMC.QT.GraphQuant",
    "category": "Type",
    "text": "GraphQuant(N::Integer, M::Integer, Γ::Float64, β::Float64, Gconstr, args...) <: DoubleGraph\n\nA DoubleGraph which implements a quantum Ising spin model in a transverse magnetic field, using the Suzuki-Trotter transformation. This allows to model the quantum transverse field case for any classical Ising model previously defined. This kind of graph can be simulated efficiently with rrrMC.\n\nN is the number of spins, M the number of Suzuki-Trotter replicas, Γ the transverse field, β the inverse temperature. GConstr is the (classical) graph constructor, and args the arguments to the contructor.\n\nSee also Qenergy and transverse_mag.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#Quantum-models-with-transverse-fields-1",
    "page": "Built-in graphs",
    "title": "Quantum models with transverse fields",
    "category": "section",
    "text": "GraphQuant"
},

{
    "location": "graphs-builtin.html#RRRMC.RE.GraphRobustEnsemble",
    "page": "Built-in graphs",
    "title": "RRRMC.RE.GraphRobustEnsemble",
    "category": "Type",
    "text": "GraphRobustEnsemble(N::Integer, M::Integer, γ::Float64, β::Float64, Gconstr, args...) <: DoubleGraph\n\nA DoubleGraph which implements a \"Robust Ensemble\" model, given any other Ising model previously defined. This kind of graph can be simulated efficiently with rrrMC.\n\nN is the number of spins, M the number of replicas, γ the interaction strength, β the inverse temperature. GConstr is the (original) graph constructor, and args the arguments to the contructor.\n\nThis is similar to GraphLocalEntropy, but the reference configuration is traced out.\n\nSee also REenergies.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#Robust-Ensemble-models-1",
    "page": "Built-in graphs",
    "title": "Robust Ensemble models",
    "category": "section",
    "text": "GraphRobustEnsemble"
},

{
    "location": "graphs-builtin.html#RRRMC.LE.GraphLocalEntropy",
    "page": "Built-in graphs",
    "title": "RRRMC.LE.GraphLocalEntropy",
    "category": "Type",
    "text": "GraphLocalEntropy(N::Integer, M::Integer, γ::Float64, β::Float64, Gconstr, args...) <: DoubleGraph\n\nA DoubleGraph which implements a \"Local Entropy\" model, given any other Ising model previously defined. This simulates a replicated system in which each replica interacts with an extra \"reference\" configuration. This kind of graph can be simulated efficiently with rrrMC.\n\nN is the number of spins, M the number of replicas, γ the interaction strength, β the inverse temperature. GConstr is the (original) graph constructor, and args the arguments to the contructor.\n\nThis is similar to GraphRobustEnsemble, but the reference is simulated explicitly.\n\nSee also LEenergies and cenergy.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#Local-Entropy-models-1",
    "page": "Built-in graphs",
    "title": "Local Entropy models",
    "category": "section",
    "text": "GraphLocalEntropy"
},

{
    "location": "graphs-builtin.html#RRRMC.Empty.GraphEmpty",
    "page": "Built-in graphs",
    "title": "RRRMC.Empty.GraphEmpty",
    "category": "Type",
    "text": "GraphEmpty(N::Integer) <: SimpleGraph{Int}\n\nA trivial SimpleGraph type with N free, non-interacting spins. (The energy is always 0.)\n\nOnly useful for testing/debugging purposes.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#RRRMC.TwoSpin.GraphTwoSpin",
    "page": "Built-in graphs",
    "title": "RRRMC.TwoSpin.GraphTwoSpin",
    "category": "Type",
    "text": "GraphTwoSpin() <: DiscrGraph{Int}\n\nA trivial DiscrGraph type with 2 spins inteacting ferromagnetically (J=1), without fields.\n\nOnly useful for testing/debugging purposes.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#RRRMC.ThreeSpin.GraphThreeSpin",
    "page": "Built-in graphs",
    "title": "RRRMC.ThreeSpin.GraphThreeSpin",
    "category": "Type",
    "text": "GraphThreeSpin() <: DiscrGraph{Int}\n\nA trivial DiscrGraph type with 3 spins, ferromagnetic interactions (J=1), no fields, and periodic boundary conditions.\n\nOnly useful for testing/debugging purposes.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#RRRMC.Fields.GraphFields",
    "page": "Built-in graphs",
    "title": "RRRMC.Fields.GraphFields",
    "category": "Type",
    "text": "GraphFields(N::Integer, LEV::Tuple = (1,)) <: DiscrGraph\n\nA simple DiscrGraph type with N non-interacting variables, each of which is subject to a local field. The fields are extracted at random from LEV, which must be a Tuple of Reals.\n\nMostly useful for testing/debugging purposes.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#RRRMC.Fields.GraphFieldsNormalDiscretized",
    "page": "Built-in graphs",
    "title": "RRRMC.Fields.GraphFieldsNormalDiscretized",
    "category": "Type",
    "text": "GraphFieldsNormalDiscretized(N::Integer, LEV::Tuple) <: DoubleGraph{Float64,GraphFields}\n\nA simple DoubleGraph type with N non-interacting variables, each of which is subject to a local field. The fields are extracted independently from a normal distribution with unit variance, and then are discretized using the values in LEV, which must be a Tuple of Reals.\n\nMostly useful for testing/debugging purposes.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#RRRMC.Ising1D.GraphIsing1D",
    "page": "Built-in graphs",
    "title": "RRRMC.Ising1D.GraphIsing1D",
    "category": "Type",
    "text": "GraphIsing1D(N::Integer) <: DiscrGraph{Int}\n\nA simple 1-dimensional DiscrGraph type with N spins, antiferromagnetic interactions (J=-1), no fields, and periodic boundary conditions.\n\nMostly useful for testing/debugging purposes.\n\n\n\n"
},

{
    "location": "graphs-builtin.html#Trivial-models-used-for-testing-and-debugging-1",
    "page": "Built-in graphs",
    "title": "Trivial models used for testing and debugging",
    "category": "section",
    "text": "GraphEmptyGraphTwoSpinGraphThreeSpinGraphFieldsGraphFieldsNormalDiscretizedGraphIsing1D"
},

{
    "location": "interface.html#",
    "page": "Graphs interface",
    "title": "Graphs interface",
    "category": "page",
    "text": ""
},

{
    "location": "interface.html#interface-1",
    "page": "Graphs interface",
    "title": "Graphs interface",
    "category": "section",
    "text": "CurrentModule = RRRMCThis page contains all the functions which are needed when implementing a graph type. See the built-in graphs for concrete examples (in particular, the SK, RRG and EA family of graphs have the most complete implementations). See also the documentation for the Config type."
},

{
    "location": "interface.html#RRRMC.Interface.energy",
    "page": "Graphs interface",
    "title": "RRRMC.Interface.energy",
    "category": "Function",
    "text": "energy(X::AbstractGraph, C::Config)\n\nReturns the energy of graph X in the configuration C. This is always invoked at the beginning of standardMC, rrrMC, bklMC and wtmMC. Subsequently, delta_energy is used instead.\n\nAll graphs must implement this function.\n\nIt must also be used to initialize/reset the cache for a given graph, if any (see update_cache!).\n\n\n\n"
},

{
    "location": "interface.html#RRRMC.Interface.delta_energy",
    "page": "Graphs interface",
    "title": "RRRMC.Interface.delta_energy",
    "category": "Function",
    "text": "delta_energy(X::AbstractGraph, C::Config, move::Int)\n\nReturns the energy difference that would be associated to flipping the spin move.\n\nA default fallback implementation based on energy is provided, to be used for debugging, but having an efficient implementation for each graph is critical for performance.\n\nNote: when X is a DiscrGraph, the absolute value of the result must be contained in the tuple returned by allΔE – no approximations are allowed, and missing values will cause crashes (unless Julia is run with the --check-bounds=yes option, in which case they will cause errors).\n\nNote: this function is always invoked before performing the flip, unlike in update_cache! and update_cache_residual!.\n\n\n\n"
},

{
    "location": "interface.html#RRRMC.Interface.update_cache!",
    "page": "Graphs interface",
    "title": "RRRMC.Interface.update_cache!",
    "category": "Function",
    "text": "update_cache!(X::AbstractGraph, C::Config, move::Int)\n\nA function which is called every time a spin is flipped. This may happen:\n\nwhen a move is accepted, in standardMC, rrrMC, bklMC and wtmMC.\nwhen a move is attempted to evaluate the effect on the neighbors, in rrrMC.\n\nmove is the spin index. By default, this function does nothing, but it may be overloaded by particular graph types.\n\nWhen X is a DoubleGraph, there is a default implementation which first calls update_cache! on inner_graph(X), then calls update_cache_residual! on X.\n\nNote: this function is always invoked after the flip has been performed, unlike in delta_energy and delta_energy_residual.\n\n\n\n"
},

{
    "location": "interface.html#RRRMC.Interface.getN",
    "page": "Graphs interface",
    "title": "RRRMC.Interface.getN",
    "category": "Function",
    "text": "getN(X::AbstractGraph)\n\nReturns the number of spins for a graph. The default implementation just returns X.N.\n\n\n\n"
},

{
    "location": "interface.html#RRRMC.Interface.neighbors",
    "page": "Graphs interface",
    "title": "RRRMC.Interface.neighbors",
    "category": "Function",
    "text": "neighbors(X::AbstractGraph, i::Int)\n\nReturns an iterable with all the neighbors of spin i. This is required by rrrMC, bklMC and wtmMC since those methods need to evaluate the effect of flipping a spin on its neighbors' local fields. It is not required by standardMC.\n\nFor performance reasons, it is best if the returned value is stack-allocated rather than heap-allocated, e.g. it is better to return an NTuple than a Vector.\n\n\n\n"
},

{
    "location": "interface.html#Functions-used-by-all-graph-types-1",
    "page": "Graphs interface",
    "title": "Functions used by all graph types",
    "category": "section",
    "text": "energydelta_energyupdate_cache!getNneighbors"
},

{
    "location": "interface.html#RRRMC.Interface.allΔE",
    "page": "Graphs interface",
    "title": "RRRMC.Interface.allΔE",
    "category": "Function",
    "text": "allΔE{P<:DiscrGraph}(::Type{P})\n\nReturns a tuple of all possible non-negative values that can be returned by delta_energy. This must be implemented by all DiscrGraph objects in order to use rrrMC or bklMC.\n\nFor performance reasons, it is best if the result can be computed from the type of the graph alone (possibly using a generated function).\n\n\n\n"
},

{
    "location": "interface.html#Functions-used-by-DiscrGraph-models-1",
    "page": "Graphs interface",
    "title": "Functions used by DiscrGraph models",
    "category": "section",
    "text": "allΔE"
},

{
    "location": "interface.html#RRRMC.Interface.inner_graph",
    "page": "Graphs interface",
    "title": "RRRMC.Interface.inner_graph",
    "category": "Function",
    "text": "inner_graph(X::DoubleGraph)\n\nReturns the internal graph used by the given DoubleGraph. The default implementation simply returns X.X0.\n\n\n\n"
},

{
    "location": "interface.html#RRRMC.Interface.delta_energy_residual",
    "page": "Graphs interface",
    "title": "RRRMC.Interface.delta_energy_residual",
    "category": "Function",
    "text": "delta_energy_residual(X::DoubleGraph, C::Config, move::Int)\n\nReturns the residual part of the energy difference produced if the spin move would be flipped, excluding the contribution from the internal SimpleGraph (see inner_graph).\n\nSee also delta_energy. There is a default fallback implementation, but it should be overloaded for efficiency.\n\n\n\n"
},

{
    "location": "interface.html#RRRMC.Interface.update_cache_residual!",
    "page": "Graphs interface",
    "title": "RRRMC.Interface.update_cache_residual!",
    "category": "Function",
    "text": "update_cache_residual!(X::DoubleGraph, C::Config, move::Int)\n\nCalled internally by the default update_cache! when the argument is a DoubleGraph. Can be useful to overload this if the residual part of the graph has an indipendent cache.\n\nBy default, it does nothing.\n\n\n\n"
},

{
    "location": "interface.html#Functions-used-by-DoubleGraph-models-1",
    "page": "Graphs interface",
    "title": "Functions used by DoubleGraph models",
    "category": "section",
    "text": "inner_graphdelta_energy_residualupdate_cache_residual!"
},

{
    "location": "interface.html#RRRMC.QT.Qenergy",
    "page": "Graphs interface",
    "title": "RRRMC.QT.Qenergy",
    "category": "Function",
    "text": "Qenergy(X::DoubleGraph, C::Config)\n\nWhen using the Suzuki-Trotter transformation to simulate quantum systems in a transverse magnetic field with a replicated classical system, this function should be used to obtain the average value of the Hamiltonian observable (divided by the number of spins).\n\n\n\n"
},

{
    "location": "interface.html#RRRMC.QT.transverse_mag",
    "page": "Graphs interface",
    "title": "RRRMC.QT.transverse_mag",
    "category": "Function",
    "text": "transverse_mag(X::DoubleGraph, C::Config, β::Float64)\n\nWhen using the Suzuki-Trotter transformation to simulate quantum systems in a transverse magnetic field with a replicated classical system, this function should be used to obtain the average value of the transverse magnetization observable.\n\n\n\n"
},

{
    "location": "interface.html#Functions-specific-to-quantum-models-1",
    "page": "Graphs interface",
    "title": "Functions specific to quantum models",
    "category": "section",
    "text": "Qenergytransverse_mag"
},

{
    "location": "interface.html#RRRMC.RE.REenergies",
    "page": "Graphs interface",
    "title": "RRRMC.RE.REenergies",
    "category": "Function",
    "text": "REenergies(X::GraphRobustEnsemble)\n\nReturns a Vector with the individual energy (as defined by the original model) of each replica in a GraphRobustEnsemble graph.\n\n\n\n"
},

{
    "location": "interface.html#Functions-specific-to-robust-ensemble-models-1",
    "page": "Graphs interface",
    "title": "Functions specific to robust-ensemble models",
    "category": "section",
    "text": "REenergies"
},

{
    "location": "interface.html#RRRMC.LE.LEenergies",
    "page": "Graphs interface",
    "title": "RRRMC.LE.LEenergies",
    "category": "Function",
    "text": "LEenergies(X::GraphLocalEntropy)\n\nReturns a Vector with the individual energy (as defined by the original model) of each replica in a GraphLocalEntropy graph.\n\n\n\n"
},

{
    "location": "interface.html#RRRMC.LE.cenergy",
    "page": "Graphs interface",
    "title": "RRRMC.LE.cenergy",
    "category": "Function",
    "text": "cenergy(X::GraphLocalEntropy)\n\nReturns a the individual energy (as defined by the original model) of the reference configuration in a GraphLocalEntropy graph.\n\n\n\n"
},

{
    "location": "interface.html#Functions-specific-to-local-entropy-models-1",
    "page": "Graphs interface",
    "title": "Functions specific to local-entropy models",
    "category": "section",
    "text": "LEenergiescenergy"
},

]}
