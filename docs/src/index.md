# RRRMC.jl documentation

```@meta
CurrentModule = RRRMC
```

This code implements the Reduced-Rejection-Rate (RRR) Monte Carlo method for Ising spin models described in the paper
*"A method to reduce the rejection rate in Monte Carlo Markov Chains"* by C. Baldassi,
J. Stat. Mech. Theor. Exp., (2017) 3, 033301 [doi:10.1088/1742-5468/aa5335](https://doi.org/10.1088/1742-5468/aa5335)
([arXiv](http://arxiv.org/abs/1608.05899)).

It also provides:
* a standard Metropolis-Hastings sampler
* a generalized implementation of the BKL method described in the paper
  ["A new algorithm for Monte Carlo simulation of Ising spin systems"](https://doi.org/10.1016/0021-9991(75)90060-1) by A.B. Bortz, M.H. Kalos and J.L. Lebowitz.
  The generalization consists in not requiring that the energy shifts are discrete.
* an implementation of the Waiting time method described in the paper
  ["Faster Monte Carlo simulations at low temperatures. The waiting time method"](https://doi.org/10.1016/S0010-4655(01)00412-X) by J. Dall and P. Sibani.
* an implementation of the "τ-Extremal Optimization" heuristic technique described in the paper
  ["Optimization with Extremal Dynamics"](https://doi.org/10.1103/PhysRevLett.86.5211) by S. Boettcher and A. G. Percus.

The code is written in [Julia](http://julialang.org). It requires Julia `1.0`.

## Installation

To install the package, use Julia's package manager: from the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
(v1.3) pkg> add RRRMC
```

Or, equivalently, via the Pkg API:

```
julia> import Pkg; Pkg.add("RRRMC")
```

Dependencies will be installed automatically.

## Usage

The module is loaded as any other Julia module:

```
julia> using RRRMC
```

The module provides four functions which implement Monte Carlo Markov Chain algorithms on Ising spin models:

* [`standardMC`](@ref): a standard Metropolis-Hastings sampler
* [`rrrMC`](@ref): the reduced-rejection-rate (RRR) method
* [`bklMC`](@ref): the Bortz-Kalos-Lebowitz (BKL) method
* [`wtmMC`](@ref): the waiting-time method (WTM)

The interface for these four algorithms is documented in the [Sampling algorithms](@ref algorithms) page, and it
is essentially identical: they take as arguments a graph, an inverse temperature parameter `β`, and the number of
Monte Carlo iterations to perform (or, for `wtmMC`, of samples to collect). However, the sampling methodology changes
based on the type of model, see the [Graph types](@ref graphtype) page.

The module also provides the function [`extremal_opt`](@ref), which implements the τ-Extremal Optimization heuristic
technique, used to seek the ground state of a model. The interface for this function is similar, but not identical,
to the others, due to its different purpose.

All of the above functions allow accessing the internal state during the iteration at regular intervals, via the `hook` keyword
argument. They also return the final configuration of the system, which is stored in an object of type
[`Config`](@ref). (The `extremal_opt` function also returns the best configuration found.)

The code comes with some [built-in graphs](@ref builtin), but provides an [interface](@ref interface) to write
user-defined models.

!!! note

    The four sampling functions are the only names exported by the module;
    all other function and types must be qualified with the `RRRMC` module
    name.

## Manual

```@contents
Pages = [
    "algorithms.md",
    "graph-types.md",
    "graphs-builtin.md",
    "interface.md",
    ]
Depth = 3
```
