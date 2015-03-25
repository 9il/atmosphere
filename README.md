[![Dub](https://img.shields.io/badge/dub-code.dlang.org-FF4081.svg)](http://code.dlang.org/packages/atmosphere_gm)

Atmosphere GM contains Maximum Likelihood Estimation algorithms, density functions, random observations generators, etc.

# Travis-CI Status
[![Coverage Status](https://coveralls.io/repos/9il/atmosphere_gm/badge.svg?branch=master)](https://coveralls.io/r/9il/atmosphere_gm?branch=master)
[![Join the chat at https://gitter.im/9il/atmosphere_gm](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/9il/atmosphere_gm?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
+ [![Build Status](https://travis-ci.org/9il/atmosphere_gm.svg)](https://travis-ci.org/9il/atmosphere_gm) - [Atmosphere GM](https://travis-ci.org/9il/atmosphere_gm)
+ [![Build Status](https://travis-ci.org/9il/atmosphere_gm_test.svg)](https://travis-ci.org/9il/atmosphere_gm_test) - [Atmosphere GM Test](https://travis-ci.org/9il/atmosphere_gm_test)

#Features
 1. Normal variance-mean mixtures
  + Generalized hyperbolic distribution
  + Generalized variance-gamma distribution
 2. Separating mixtures of probability distributions
  + Grid methods
  + Likelihood maximization
  + Optimization over sliding window
 3. Generalized probability distributions
  + Density functions
  + Cumulative functions
  + Quantiles
  + Random observations generators
  + Maximum Likelihood Estimations (MLE)

<img src="http://9il.github.io/atmosphere_gm/doc/images/GHyp_0148.svg" alt="Generalized hyperbolic distribution" width="280" />
<img src="http://9il.github.io/atmosphere_gm/doc/images/GV-gamma_0120.svg" alt="Generalized variance-gamma distribution" width="280" />

#### What can I use this package for? 
Atmosphere GM can be used for risk management in economics, finance and thermonuclear reactors ;-)

#Documentation
Documentation (API) can be found [here](http://9il.github.io/atmosphere_gm/doc/atmosphere.html).

#Installation
## BLAS & LAPACK
You need BLAS and LAPACL libraries installed.

If you're on **Ubuntu**, you can install default package 
```
sudo apt-get install libblas-dev liblapack-dev 
```

**OS X** comes with the [Accelerate framework](https://developer.apple.com/library/mac/documentation/Accelerate/Reference/BLAS_Ref/index.html#//apple_ref/doc/uid/TP40009457) built in. 

There is no instruction for **Windows**. You are welcome to create it.

## Intro to D
1. Install D [compiler](http://dlang.org/download.html). ~~The [2.067-b2](https://dlang.dawg.eu/downloads/dmd.2.067.0-b2/) or later is required~~.
2. Install [DUB registry](http://code.dlang.org/download)
3. Read about [DUB](http://code.dlang.org/about)
4. Read about [DUB package format](http://code.dlang.org/package-format)

## Atmosphere GM
To use this package put [the dependency]((http://code.dlang.org/packages/atmosphere_gm)) into your project's
[dub](http://code.dlang.org/about).json into the dependencies section
and the following imports into your program
```D
import atmosphere;
```
If you want to write library use detailed imports
```D
import atmosphere.pdf;
import atmosphere.estimate.generalized_inverse_gaussian;
import atmosphere.finitemixture;
import atmosphere.mixture : MixtureOptimizer, MixtureOptimizerException;
```

# Compilers and optimization
The [DMD](http://dlang.org/download.html) compiler is easy way to start.
To compile your program in release mode use the following build options
```shell
dub build --build=release
```
## LDC
It is suggested the [LLVM D Compiler](https://github.com/ldc-developers/ldc/releases) be used for benchmarks.
To compile your program in release mode with LDC use the following build options
```
dub build --build=release --compiler=ldc2
```
To fine-tune your program for native CPU add the following code into your `dub.json`:
```
{
	...
	"dflags-ldc": ["-mcpu=native"],
}
```
For more options run `ldc2 -help`.

See also [Atmosphere GM Test](https://github.com/9il/atmosphere_gm_test). 

# TODO & Contribution
Contribution is welcome.
#### TODO list
+ **Windows** instruction for BLAS & LAPACK installation.
+ Parameter estimation algorithms for generalized inverse Gaussian and generalized gamma samples.
+ Publication references
+ More `unittest`s
+ Examples
+ [ReadTheDocs](https://readthedocs.org) documentation.
