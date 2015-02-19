/++
This module publicly imports all classes.

Import:
--------
import atmosphere;
--------

Class_hierarchy:
- [MixtureOptimizer](atmosphere/mixture/MixtureOptimizer.html)
	- [ExpectationMaximization](atmosphere/mixture/ExpectationMaximization.html)
		- [LikelihoodAscentEM](atmosphere/mixture/LikelihoodAscentEM.html) : [LikelihoodAscent](atmosphere/mixture/LikelihoodAscent.html)
	- [GradientDescent](atmosphere/mixture/GradientDescent.html)
	- [GradientDescentPartial](atmosphere/mixture/GradientDescentPartial.html)
		- [LikelihoodAscentGradient](atmosphere/mixture/LikelihoodAscentGradient.html) : [LikelihoodAscent](atmosphere/mixture/LikelihoodAscent.html)
	- [CoordinateDescent](atmosphere/mixture/CoordinateDescent.html)
	- [CoordinateDescentPartial](atmosphere/mixture/CoordinateDescentPartial.html)
		- [LikelihoodAscentCoordinate](atmosphere/mixture/LikelihoodAscentCoordinate.html) : [LikelihoodAscent](atmosphere/mixture/LikelihoodAscent.html)
	- [NvmmLikelihoodAscentEM](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEM.html) : [LikelihoodAscent](atmosphere/mixture/LikelihoodAscent.html)
		- [NvmmLikelihoodAscentEMEM](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEMEM.html)
		- [NvmmLikelihoodAscentEMGradient](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEMGradient.html)
		- [NvmmLikelihoodAscentEMCoordinate](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEMCoordinate.html)

See_Also: [Atmosphere GM Test](https://github.com/9il/atmosphere_gm_test)

+/

/+

Table of separating mixtures algorithms:
	$(TABLE
		$(TR  $(TH Class) $(TH BLAS Level 2 % of all computations) $(TH Non-BLAS parallel) $(TH   Speed (one thread & big data) ) $(TH Global maximum for convex function) $(TH Parameterized))
		$(TR $(TD [GradientDescent](atmosphere/mixture/GradientDescent.html), number of mixture components > ~16) $(TD ~95-99%) $(TD No) $(TD  normal ) $(TD Yes) $(TD No) )
		$(TR $(TD [GradientDescent](atmosphere/mixture/GradientDescent.html), number of mixture components < ~16) $(TD ~95-99%) $(TD No) $(TD  fast ) $(TD Yes) $(TD No) )
		$(TR $(TD [CoordinateDescent](atmosphere/mixture/CoordinateDescent.html)) $(TD ~10-50%) $(TD No) $(TD  fast ) $(TD Yes) $(TD No) )
		$(TR $(TD [CoordinateDescentPartial](atmosphere/mixture/LikelihoodAscentCoordinate.html)) $(TD ~10-50%) $(TD No) $(TD  fast, always faster then CoordinateDescent ) $(TD Yes) $(TD No) )
		$(TR $(TD [NvmmLikelihoodAscentEMEM](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEMEM.html)) $(TD ~10%) $(TD Yes) $(TD  slow ) $(TD No) $(TD Yes) )
		$(TR $(TD [NvmmLikelihoodAscentEMGradient](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEMGradient.html)) $(TD ~10%) $(TD Yes) $(TD  slow ) $(TD No) $(TD Yes) )
		$(TR $(TD [NvmmLikelihoodAscentEMCoordinate](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEMCoordinate.html)) $(TD ~5%) $(TD Partial) $(TD  slow ) $(TD No) $(TD Yes) ) )
+/
/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere;

///
unittest
{
	import core.time;
	import std.random;
	import std.range;
	import std.stdio;
	import atmosphere;

	alias F = double;

	final class ProperGeneralizedInverseGaussianQuantile(T) : NumericQuantile!T {
		this(T lambda, T eta, T omega) {
			auto cdf = new ProperGeneralizedInverseGaussianCDF!T(lambda, eta, omega);
			super(cdf, -1000, 1000);	
		}
	}

	nothrow @nogc bool findRootTolerance(F a, F b) { return b/a < 1.001;}
	
	immutable quantileL  = 0.01;
	immutable quantileR  = 0.99;
	immutable gridSize   = 100;
	immutable dur        = TickDuration.from!"msecs"(100);
	immutable lambda     = 2;
	immutable eta        = 1; 
	immutable omega      = 2.3;
	immutable beta       = 0.5;
	immutable sampleSize = 1000;
	// GIG quantile function
	auto qf              = new ProperGeneralizedInverseGaussianQuantile!F(lambda, eta, omega);
	// GHyp random number generator
	auto rng             = new ProperGeneralizedHyperbolicRNG!F(rndGen, lambda, eta, omega, beta);
	// left GIG bound
	immutable begin      = qf(quantileL);
	// right GIG bound
	immutable end        = qf(quantileR);
	// grid's step
	immutable step       = (end-begin)/gridSize;
	// GIG grid
	immutable grid       = iota(begin, end+step/2, step).array;
	// Normal PDFs for common algorithms
	auto pdfs            = grid.map!(u => NvmmLikelihoodAscentEM!F.CorePDF(beta, u)).array;
	// GHyp sample
	immutable sample     = cast(immutable) rng.take(sampleSize).array;
	// Mixture optimizer
	auto optimizer       = new LikelihoodAscentEM!F(pdfs.length, sample.length);
	// Puts sample
	optimizer.put(pdfs, sample);
	// Tuple of time and iterations count
	immutable result     = optimizer.evaluate(dur, &findRootTolerance);
}

public import atmosphere.cdf;
public import atmosphere.derivative;
public import atmosphere.finitemixture;
public import atmosphere.mixture;
public import atmosphere.moment;
public import atmosphere.params;
public import atmosphere.pdf;
public import atmosphere.quantile;
public import atmosphere.random;

public import atmosphere.estimate;
public import atmosphere.likelihood;
