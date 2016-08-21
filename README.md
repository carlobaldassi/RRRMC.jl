# RRRMC.jl

[![Build Status](https://travis-ci.org/carlobaldassi/RRRMC.jl.svg?branch=master)](https://travis-ci.org/carlobaldassi/RRRMC.jl)

[![Coverage Status](https://coveralls.io/repos/carlobaldassi/RRRMC.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/carlobaldassi/RRRMC.jl?branch=master)

[![codecov.io](http://codecov.io/github/carlobaldassi/RRRMC.jl/coverage.svg?branch=master)](http://codecov.io/github/carlobaldassi/RRRMC.jl?branch=master)

This code implements the Reduced-Rejection-Rate (RRR) Monte Carlo method for Ising spin models described in the paper
"A method to reduce the rejection rate in Monte Carlo Markov Chains on Ising spin models" by C. Baldassi.

It also provides a standard Metropolis-Hastings sampler, and an implementation of the BKL method described in the paper
["A new algorithm for Monte Carlo simulation of Ising spin systems"][BKLpaper] by A.B. Bortz, M.H. Kalos and J.L. Lebowitz
 	
The code is written in [Julia].

The package was tested on Julia `0.4` and `0.5`.

### Installation

To install the module, use this command from within Julia:

```
julia> Pkg.clone("https://github.com/carlobaldassi/RRRMC.jl")
```

### Documentation

The documentation is currently under construction. At the moment, the functions documentation is available at the Julia REPL,
after the module is installed and loaded: press `?` to enter the Julia help mode, then look at the documentation for the functions
`standardMC`, `bklMC` and `rrrMC`, and the other functions and types mentioned therein.

[Julia]: http://julialang.org
[BKLpaper]: http://www.sciencedirect.com/science/article/pii/0021999175900601
 	
 	
 	
 	
 	
 	
 	
 	
 	
 	
 	
