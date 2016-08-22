# RRRMC.jl documentation

```@meta
CurrentModule = RRRMC
```

This package implements the Reduced-Rejection-Rate (RRR) Monte Carlo method for Ising spin models described in the paper
"A method to reduce the rejection rate in Monte Carlo Markov Chains on Ising spin models" by C. Baldassi.

It also provides a standard Metropolis-Hastings sampler, and an implementation of the BKL method described in the paper
["A new algorithm for Monte Carlo simulation of Ising spin systems"](http://www.sciencedirect.com/science/article/pii/0021999175900601) by A.B. Bortz, M.H. Kalos and J.L. Lebowitz

The code is written in [Julia](http://julialang.org), and tested against Julia `0.4`, `0.5` and *current* `0.6-dev` on Linux, OS X, and Windows.

### Installation

To install the module, use this command from within Julia:

```
julia> Pkg.clone("https://github.com/carlobaldassi/RRRMC.jl")
```

Dependencies will be installed automatically.

### Usage

The module is loaded as any other Julia module:

```
julia> using RRRMC
```

The module provides three functions which implement Monte Carlo Markov Chain algorithms on Ising spin models:

* [`standardMC`](@ref): a standard Metropolis-Hastings sampler
* [`rrrMC`](@ref): the reduced-rejection-rate (RRR) method
* [`bklMC`](@ref): the Bortz-Kalos-Lebowitz (BKL) method

The interface for these three algorithms is documented in the [Sampling algorithms](@ref algorithms) page, and it
is essentially identical: they take as arguments a graph, an inverse temperature parameter `β`, and the number of
Monte Carlo iterations to perform. However, `rrrMC` and `bklMC` can only be used on some type of models, see the
[Graph types](@ref graphtype) page.

These functions allow accessing the internal state during the iteration at regular intervals, via the `hook` keyword
argument. They also return the final configuration of the system, which is stored in an object of type
[`Config`](@ref).

The code comes with some [built-in graphs](@ref builtin), but provides an [interface](@ref interface) to write
user-defined models.

!!! note
    The three sampling functions are the only names exported by the module;
    all other function and types must be qualified with the `RRRMC` module
    name.

### Manual

```@contents
Pages = [
    "algorithms.md",
    "graph-types.md",
    "graphs-builtin.md",
    "interface.md",
    ]
Depth = 3
```
