# RRRMC.jl

| **Documentation**                       | **Build Status**                                                                                |
|:---------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-latest-img]][docs-latest-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] |

This code implements the Reduced-Rejection-Rate (RRR) Monte Carlo method for Ising spin models described in the paper
"A method to reduce the rejection rate in Monte Carlo Markov Chains on Ising spin models" by C. Baldassi.

It also provides a standard Metropolis-Hastings sampler, and an implementation of the BKL method described in the paper
["A new algorithm for Monte Carlo simulation of Ising spin systems"][BKLpaper] by A.B. Bortz, M.H. Kalos and J.L. Lebowitz

The code is written in [Julia], and tested against Julia `0.4`, `0.5` and *current* `0.6-dev` on Linux, OS X, and Windows.

### Installation

To install the module, use this command from within Julia:

```
julia> Pkg.clone("https://github.com/carlobaldassi/RRRMC.jl")
```

### Documentation

- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

[Julia]: http://julialang.org
[BKLpaper]: http://www.sciencedirect.com/science/article/pii/0021999175900601

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://carlobaldassi.github.io/RRRMC.jl/latest

[travis-img]: https://travis-ci.org/carlobaldassi/RRRMC.jl.svg?branch=master
[travis-url]: https://travis-ci.org/carlobaldassi/RRRMC.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/bq8jj4u0dx6x6xm1/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/carlobaldassi/rrrmc-jl/branch/master

[codecov-img]: https://codecov.io/gh/carlobaldassi/RRRMC.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/carlobaldassi/RRRMC.jl
