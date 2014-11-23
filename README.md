#Features
  0. Separating mixtures of probability distributions
  1. Grid methods
  2. Parameterized algorithms
  3. Optimization over sliding window

#Documentation
Documentation can be found [here](http://9il.github.io/atmosphere_gm/atmosphere.package.html).

#Installation
To use [this package](http://code.dlang.org/packages/atmosphere_gm), put the following dependency into your project's
[dub](http://code.dlang.org/about).json into the dependencies section:
```json
{
	...
	"dependencies": {
		"atmosphere_gm": ">=0.0.2"
	}
}
```

##BLAS
###ubuntu
```shell
sudo apt-get install libblas-dev
```
### OS X
Mac OS X comes with the Accelerate framework built in.

### Windows
TODO

##Intro to D
1. Install D [compiler](http://dlang.org/download.html)
2. Install [DUB registry](http://code.dlang.org/download)
3. Read about [DUB](http://code.dlang.org/about)
4. Read about [DUB package format](http://code.dlang.org/package-format)
5. Start with [example](https://github.com/9il/atmosphere_gm/tree/master/examples/normal_variance_mean_mixture)

#Benchmarking
It is suggested the [llvm D compiler](https://github.com/ldc-developers/ldc/releases) be used for benchmarks.
Project requires LDC version >= 0.15.0 or DMD >= 2.066, or corresponding GDC release.
The [DMD](http://dlang.org/download.html) is easy way to start.

#Changelog
##v0.0.2
1. Exceptions added; #3
2. NVMM test updated to show `grid`; #6
