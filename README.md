# RRRMC.jl

| **Documentation**                       | **PackageEvaluator**                                    | **Build Status**                                                                                | **Releases**                     |
|:---------------------------------------:|:-------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:--------------------------------:|
| [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.5-img]][pkg-url] [![][pkg-0.6-img]][pkg-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] | [![DOI][zenodo-img]][zenodo-url] |

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

The code is written in [Julia], and tested against Julia `0.5`, `0.6` and *current* (at the time of writing) `0.7-DEV` on
Linux, OS X, and Windows.

### Installation

To install the module, use Julia's package manager:

```
julia> Pkg.add("RRRMC")
```

Dependencies will be installed automatically.

### Documentation

- [**STABLE**][docs-stable-url] &mdash; stable version of the documentation
- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

[Julia]: https://julialang.org
[RRRpaper]: https://doi.org/10.1088/1742-5468/aa5335
[RRRarXiv]: http://arxiv.org/abs/1608.05899
[BKLpaper]: https://doi.org/10.1016/0021-9991(75)90060-1
[WTMpaper]: https://doi.org/10.1016/S0010-4655(01)00412-X

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://carlobaldassi.github.io/RRRMC.jl/stable
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://carlobaldassi.github.io/RRRMC.jl/latest

[travis-img]: https://travis-ci.org/carlobaldassi/RRRMC.jl.svg?branch=master
[travis-url]: https://travis-ci.org/carlobaldassi/RRRMC.jl

[pkg-0.5-img]: http://pkg.julialang.org/badges/RRRMC_0.5.svg
[pkg-0.6-img]: http://pkg.julialang.org/badges/RRRMC_0.6.svg
[pkg-url]: http://pkg.julialang.org/?pkg=RRRMC

[appveyor-img]: https://ci.appveyor.com/api/projects/status/bq8jj4u0dx6x6xm1/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/carlobaldassi/rrrmc-jl/branch/master

[codecov-img]: https://codecov.io/gh/carlobaldassi/RRRMC.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/carlobaldassi/RRRMC.jl

[zenodo-img]: https://zenodo.org/badge/66179142.svg
[zenodo-url]: https://zenodo.org/badge/latestdoi/66179142
