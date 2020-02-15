# RRRMC.jl

| **Documentation**                 | **Build Status**                                              | **Releases**                     |
|:---------------------------------:|:-------------------------------------------------------------:|:--------------------------------:|
| [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url] [![][codecov-img]][codecov-url] | [![DOI][zenodo-img]][zenodo-url] |

This code implements the Reduced-Rejection-Rate (RRR) Monte Carlo method for Ising spin models described in the paper
*"A method to reduce the rejection rate in Monte Carlo Markov Chains"* by C. Baldassi,
J. Stat. Mech. Theor. Exp., (2017) 3, 033301 [doi:10.1088/1742-5468/aa5335][RRRpaper] ([arXiv][RRRarXiv]).

It also provides:
* a standard Metropolis-Hastings sampler
* a generalized implementation of the BKL method described in the paper
  ["A new algorithm for Monte Carlo simulation of Ising spin systems"][BKLpaper] by A.B. Bortz, M.H. Kalos and J.L. Lebowitz.
  The generalization consists in not requiring that the energy shifts are discrete.
* an implementation of the Waiting time method described in the paper
  ["Faster Monte Carlo simulations at low temperatures. The waiting time method"][WTMpaper] by J. Dall and P. Sibani.
* an implementation of the "τ-Extremal Optimization" heuristic technique described in the paper
  ["Optimization with Extremal Dynamics"][EOpaper] by S. Boettcher and A. G. Percus.

The code is written in [Julia]. It requires Julia `1.0`.

### Installation

To install the package, use Julia's package manager: from the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
(v1.3) pkg> add RRRMC
```

Or, equivalently, via the Pkg API:

```
julia> import Pkg; Pkg.add("RRRMC")
```

Dependencies will be installed automatically.

### Documentation

- [**STABLE**][docs-stable-url] &mdash; stable version of the documentation
- [**DEV**][docs-dev-url] &mdash; *in-development version of the documentation.*

[Julia]: https://julialang.org
[RRRpaper]: https://doi.org/10.1088/1742-5468/aa5335
[RRRarXiv]: http://arxiv.org/abs/1608.05899
[BKLpaper]: https://doi.org/10.1016/0021-9991(75)90060-1
[WTMpaper]: https://doi.org/10.1016/S0010-4655(01)00412-X
[EOpaper]: https://doi.org/10.1103/PhysRevLett.86.5211

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://carlobaldassi.github.io/RRRMC.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://carlobaldassi.github.io/RRRMC.jl/dev

[travis-img]: https://travis-ci.org/carlobaldassi/RRRMC.jl.svg?branch=master
[travis-url]: https://travis-ci.org/carlobaldassi/RRRMC.jl

[codecov-img]: https://codecov.io/gh/carlobaldassi/RRRMC.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/carlobaldassi/RRRMC.jl

[zenodo-img]: https://zenodo.org/badge/66179142.svg
[zenodo-url]: https://zenodo.org/badge/latestdoi/66179142
