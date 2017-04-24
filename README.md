# RRRMC.jl

| **Documentation**                       | **Build Status**                                                                                | **Releases**                     |
|:---------------------------------------:|:-----------------------------------------------------------------------------------------------:|:--------------------------------:|
| [![][docs-latest-img]][docs-latest-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] | [![DOI][zenodo-img]][zenodo-url] |

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

The code is written in [Julia], and tested against Julia `0.5` and *current* `0.6-pre` on Linux, OS X, and Windows.

### Installation

To install the module, use this command from within Julia:

```
julia> Pkg.clone("https://github.com/carlobaldassi/RRRMC.jl")
```

Dependencies will be installed automatically.

### Documentation

- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

[Julia]: https://julialang.org
[RRRpaper]: https://doi.org/10.1088/1742-5468/aa5335
[RRRarXiv]: http://arxiv.org/abs/1608.05899
[BKLpaper]: http://www.sciencedirect.com/science/article/pii/0021999175900601
[WTMpaper]: http://www.sciencedirect.com/science/article/pii/S001046550100412X

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://carlobaldassi.github.io/RRRMC.jl/latest

[travis-img]: https://travis-ci.org/carlobaldassi/RRRMC.jl.svg?branch=master
[travis-url]: https://travis-ci.org/carlobaldassi/RRRMC.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/bq8jj4u0dx6x6xm1/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/carlobaldassi/rrrmc-jl/branch/master

[codecov-img]: https://codecov.io/gh/carlobaldassi/RRRMC.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/carlobaldassi/RRRMC.jl

[zenodo-img]: https://zenodo.org/badge/66179142.svg
[zenodo-url]: https://zenodo.org/badge/latestdoi/66179142
