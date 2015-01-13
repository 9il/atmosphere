Atmosphere GM contains an experimental algorithms for new generation of likelihood maximization methods. The goal of this library is to create an reliable and simple in use technology that can be used by software engineers. Package is growing step by step with mathematical publications.

# Travis-CI Status
+ [![Build Status](https://travis-ci.org/9il/atmosphere_gm.svg)](https://travis-ci.org/9il/atmosphere_gm) - [Atmosphere GM](https://travis-ci.org/9il/atmosphere_gm)
+ [![Build Status](https://travis-ci.org/9il/atmosphere_gm_test.svg)](https://travis-ci.org/9il/atmosphere_gm_test) - [Atmosphere GM Test](https://travis-ci.org/9il/atmosphere_gm_test)

#Features
 1. Normal variance-mean mixtures
  + Generalized hyperbolic distribution
  + Generalized variance-gamma distribution (in progress)
 2. Separating mixtures of probability distributions
  + Grid methods
  + Likelihood maximization
  + Optimization over sliding window
 3. Generalized probability distributions
  + Density functions
  + Cumulative functions
  + Quantiles
  + Random observations generators


#Documentation
Documentation (API) can be found [here](http://9il.github.io/atmosphere_gm/doc/atmosphere.html).
## Changelog
See [release notes](https://github.com/9il/atmosphere_gm/releases).

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
1. Install D [compiler](http://dlang.org/download.html)
2. Install [DUB registry](http://code.dlang.org/download)
3. Read about [DUB](http://code.dlang.org/about)
4. Read about [DUB package format](http://code.dlang.org/package-format)

## Atmosphere GM
To use [this package](http://code.dlang.org/packages/atmosphere_gm), put the following dependency into your project's
[dub](http://code.dlang.org/about).json into the dependencies section:
```json
{
	...
	"dependencies": {
		"atmosphere_gm": ">=0.0.10"
	}
}
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

# Bugs
See [bug list](https://github.com/9il/atmosphere_gm/labels/bug).

# TODO & Contribution
Contribution is welcome.
#### TODO list
+ **Windows** instruction for BLAS & LAPACK installation.
+ Generalized variance-gamma probability density function.
+ Parameter estimation algorithms for generalized inverse Gaussian and variance-gamma samples.
+ Publication references
+ More `unittest`s
+ Examples
+ [ReadTheDocs](https://readthedocs.org) documentation.
